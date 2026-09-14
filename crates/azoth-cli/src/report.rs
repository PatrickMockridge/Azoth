//! Rendering a result for a terminal.
//!
//! Two rules the output follows, both of which matter more than they look for a
//! tool whose whole purpose is to be believed:
//!
//! * **Assumptions are printed, not implied.** The roughness default and the
//!   velocity the CLI derived from the flow are shown every time, because the
//!   user did not supply them and cannot see them otherwise.
//! * **Warnings are printed, never summarised away.** A count would let someone
//!   read the pressure drop and skip the caveats, which is the failure mode this
//!   library is organised against.

use azoth_core::{Result, Warning};
use azoth_hydraulics::known_fittings;

use crate::pipe::{FlowUnit, PipeResult};

/// Width of the label column.
const LABEL: usize = 20;

/// Print the pipe report.
pub fn print_pipe(result: &PipeResult, flow_unit: FlowUnit, flow: f64, roughness_supplied: bool) {
    print!(
        "{}",
        render_pipe(result, flow_unit, flow, roughness_supplied)
    );
}

/// The pipe report, as text.
///
/// Separate from [`print_pipe`] so the report can be asserted on rather than read
/// off a terminal. What the output *says* is the whole of this module's job - the
/// assumptions it prints and the warnings it refuses to summarise - and none of
/// that was checkable while the only way to see it was stdout.
pub fn render_pipe(
    result: &PipeResult,
    flow_unit: FlowUnit,
    flow: f64,
    roughness_supplied: bool,
) -> String {
    let mut out = String::new();

    out.push_str("\nazoth pipe\n\n");

    section(&mut out, "inputs");
    row(
        &mut out,
        "fluid",
        &format!("{} at {} C", result.fluid, result.temperature_c),
    );
    row(
        &mut out,
        "density",
        &format!("{:.4} kg/m3", result.density.value),
    );
    row(
        &mut out,
        "viscosity",
        &format!("{:.6e} Pa*s", result.viscosity.value),
    );
    row(
        &mut out,
        "diameter",
        &format!("{} m", result.diameter.value),
    );
    row(&mut out, "length", &format!("{} m", result.length.value));
    row(
        &mut out,
        "roughness",
        &format!(
            "{:.3e} m{}",
            result.roughness_m,
            if roughness_supplied {
                ""
            } else {
                "  (default - override with --roughness)"
            }
        ),
    );
    row(
        &mut out,
        "flow",
        &format!(
            "{flow} {} = {:.6} m3/s",
            flow_unit.as_str(),
            result.flow_m3_s
        ),
    );
    row(
        &mut out,
        "velocity",
        &format!(
            "{:.6} m/s  (derived from flow and bore)",
            result.velocity.value
        ),
    );

    section(&mut out, "flow");
    row(
        &mut out,
        "reynolds number",
        &format!("{:.1}  ({})", result.re, result.regime),
    );
    row(
        &mut out,
        "relative roughness",
        &format!("{:.6e}", result.relative_roughness),
    );
    let iterations = match result.iterations {
        Some(n) => format!(", {n} iterations"),
        None => " (explicit, no iteration)".to_string(),
    };
    row(
        &mut out,
        "friction factor",
        &format!(
            "{:.8}  ({}{})",
            result.friction_factor,
            result.method.as_str(),
            iterations
        ),
    );

    section(&mut out, "pressure drop");
    if result.fittings.is_empty() {
        row(&mut out, "fittings", "none");
    } else {
        row(
            &mut out,
            "fittings",
            &format!(
                "{}  ->  K = {:.6}",
                result.fittings.join(", "),
                result.k_total
            ),
        );
    }
    row(
        &mut out,
        "straight pipe",
        &format!("{:>12.4} Pa", result.dp_straight),
    );
    if !result.fittings.is_empty() {
        row(
            &mut out,
            "fittings",
            &format!(
                "{:>12.4} Pa   ({:.1}% of total)",
                result.dp_fittings,
                100.0 * result.dp_fittings / result.dp_total
            ),
        );
    }
    row(
        &mut out,
        "total",
        &format!(
            "{:>12.4} Pa   ({:.6} kPa, {:.6} bar)",
            result.dp_total,
            result.dp_total / 1000.0,
            result.dp_total / 1e5
        ),
    );

    out.push_str(&render_warnings(&result.warnings));
    out
}

fn section(out: &mut String, title: &str) {
    out.push_str(&format!("  {title}\n"));
}

fn row(out: &mut String, label: &str, value: &str) {
    out.push_str(&format!("    {label:<LABEL$} {value}\n"));
}

/// Print the warnings, in full.
///
/// The whole list, not a count and not a summary. A caveat the user did not read
/// is a caveat that was not given.
pub fn print_warnings(warnings: &[Warning]) {
    print!("{}", render_warnings(warnings));
}

/// The warnings, as text. See [`print_warnings`].
pub fn render_warnings(warnings: &[Warning]) -> String {
    let mut out = String::new();
    if warnings.is_empty() {
        out.push_str("\n  no warnings: every range check passed\n\n");
        return out;
    }

    out.push_str(&format!("\n  {} warning(s)\n", warnings.len()));
    for warning in warnings {
        out.push_str(&format!("    [{}] {}\n", warning.code, warning.message));
        if let Some(field) = &warning.field {
            out.push_str(&format!("      about: {field}\n"));
        }
    }
    out.push_str(
        "\n  A value outside its validated range is still a value, but it has not been\n  \
         checked the way an in-range one has. Read the above before using this number.\n\n",
    );
    out
}

/// Print the known fitting ids.
///
/// # Errors
/// Propagates a malformed embedded registry.
pub fn list_fittings() -> Result<()> {
    print!("{}", render_fittings()?);
    Ok(())
}

/// The fitting registry, as text. See [`list_fittings`].
///
/// # Errors
/// Propagates a malformed embedded registry.
pub fn render_fittings() -> Result<String> {
    let ids = known_fittings()?;
    let mut out = String::new();
    out.push_str(&format!("\nknown fitting ids ({}):\n\n", ids.len()));
    for fitting in azoth_hydraulics::fittings::registry()? {
        let marker = if fitting.is_estimated() {
            "  [ESTIMATED DUMMY - not engineering data]"
        } else {
            ""
        };
        out.push_str(&format!(
            "  {:<26} {:<28} n_ld = {:<6}{marker}\n",
            fitting.id, fitting.name, fitting.n_ld
        ));
    }
    out.push_str(
        "\nEvery coefficient above is a placeholder for software testing. See\n\
         data/fittings/crane_k_factors.csv. Do not size equipment with them.\n\n",
    );
    Ok(out)
}
