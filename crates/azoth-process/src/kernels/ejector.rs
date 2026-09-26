//! `unit_ops.ejector` - a motive stream, a suction stream, and the diffuser that follows.

use azoth_core::units::{
    Pressure, Velocity, joules_per_mole, joules_per_mole_kelvin, meters_per_second, pascals,
};
use azoth_core::{AzothError, Result};
use azoth_eos::{ph_flash, ps_flash};

use crate::stream::Stream;

/// The numbers the ejector's own route reached, in SI.
///
/// **None of them is an input.** `setMixingPressure` exists and the palette does not declare it,
/// so the mixing pressure is the class's own estimate - and the four velocities are what the two
/// nozzle efficiencies, the mixing efficiency and the diffuser's estimate produced from it. A
/// reader who cannot see them is looking at an outlet state with no way to tell an ejector that
/// barely drew from one that was near its limit.
pub struct EjectorNumbers {
    /// The pressure the two streams meet at, Pa: `estimateDefaultMixingPressure`, clamped to the
    /// suction pressure.
    pub mixing_pressure: Pressure,
    /// The motive nozzle's exit velocity, m/s: `sqrt(2 dh)` at its own efficiency.
    pub motive_nozzle_velocity: Velocity,
    /// The suction nozzle's, m/s: the larger of its own `sqrt(2 dh)` and the class's blended
    /// estimate.
    pub suction_nozzle_velocity: Velocity,
    /// The mixed stream's, m/s: the momentum balance, scaled by the mixing efficiency.
    pub mixing_velocity: Velocity,
    /// The diffuser's design velocity, m/s - the class's own estimate, which takes one more
    /// `v^2/2` off before the discharge flash.
    pub diffuser_velocity: Velocity,
}

/// What an ejector did: the discharge record, and the five numbers its route reached.
pub struct EjectorOutcome {
    /// The discharge record.
    pub outlet: Stream,
    /// The mixing pressure and the four velocities, or `None` where nothing flowed - the class
    /// returns the motive stream unchanged there rather than dividing by a total of zero, and a
    /// machine that ran on no fluid reached no pressure and no velocity at all.
    pub numbers: Option<EjectorNumbers>,
}

/// The ejector's four parameters, as `unit_ops.ejector` declares them.
///
/// **The palette entry was missing one of the five the class reads.** `run` expands the
/// motive stream with `efficiencyIsentropic`, a field of its own that defaults to `0.75` and
/// has its own setter - so a `motive_nozzle_efficiency` the entry did not declare was moving
/// the answer, and the three efficiencies it did declare are the suction nozzle's, the mixing
/// chamber's and the diffuser's.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EjectorSetup {
    /// The pressure the diffuser discharges at, Pa.
    pub discharge_pressure: Pressure,
    /// The motive nozzle's isentropic efficiency.
    pub motive_nozzle_efficiency: f64,
    /// The suction nozzle's isentropic efficiency.
    pub suction_nozzle_efficiency: f64,
    /// The mixing chamber's momentum-transfer efficiency.
    pub mixing_efficiency: f64,
    /// The diffuser's pressure-recovery efficiency.
    pub diffuser_efficiency: f64,
}

/// **The bar NeqSim's own estimates are written in**, which is not this library's unit.
///
/// Three of `run`'s helpers take pressures in bara and multiply the difference by `1e5` where
/// they need pascals, and one of them floors at `0.01` - a *bar*, not a pascal. A port that
/// worked in pascals throughout would be off by five orders of magnitude on that floor, so
/// the two estimates are written in the unit the class wrote them in and converted at the
/// boundary.
const BAR: f64 = 1.0e5;

/// Expand a motive stream, entrain it with a suction stream, and diffuse the mixture.
///
/// The route is `Ejector.run`'s: each stream is flashed isentropically to the **mixing
/// pressure** the class estimates, dropped to the enthalpy its own nozzle efficiency leaves,
/// and given a velocity `sqrt(2 dh)`; the two momenta mix; the mixture is flashed to the
/// static enthalpy the mixing leaves and again across the diffuser's pressure rise; and the
/// diffuser's own design velocity takes one more `v^2/2` off before the final flash.
///
/// **The two design velocities are the class's own empirical estimates, and they decide the
/// answer.** `run` takes the larger of the enthalpy's velocity and a blended estimate
/// `0.6 sqrt(2·0.02Δp/ρ) + 0.4 (50 + 30 log10(1 + 5 Q))` clamped to `[25, 120]` m/s at the
/// suction, and a mirror pair `[10, 60]` at the diffuser outlet - so a port that skipped them
/// would be solving a different machine.
///
/// # Errors
/// [`AzothError::InvalidInput`] if an efficiency is outside `(0, 1]` or the discharge pressure
/// is below the suction's, and whatever the flashes refuse.
pub fn ejector(motive: &Stream, suction: &Stream, setup: EjectorSetup) -> Result<EjectorOutcome> {
    for (name, value) in [
        ("motive_nozzle_efficiency", setup.motive_nozzle_efficiency),
        ("suction_nozzle_efficiency", setup.suction_nozzle_efficiency),
        ("mixing_efficiency", setup.mixing_efficiency),
        ("diffuser_efficiency", setup.diffuser_efficiency),
    ] {
        if !(value > 0.0 && value <= 1.0) {
            return Err(AzothError::invalid_input(
                name,
                format!("an efficiency is a fraction in (0, 1], and {value} is not"),
            ));
        }
    }
    if setup.discharge_pressure.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "discharge_pressure",
            format!(
                "a discharge pressure is positive, and {} Pa is not",
                setup.discharge_pressure.value
            ),
        ));
    }

    let motive_mass = motive.mass_flow()?;
    let suction_mass = suction.mass_flow()?;
    let total_mass = motive_mass + suction_mass;
    let suction_pressure_bar = suction.p.value / BAR;
    let discharge_pressure_bar = setup.discharge_pressure.value / BAR;

    // `run` returns the motive stream unchanged when nothing flows, rather than dividing by
    // the total.
    if total_mass <= 0.0 {
        return Ok(EjectorOutcome {
            outlet: motive.clone(),
            numbers: None,
        });
    }

    let mut mixing_pressure_bar = estimate_mixing_pressure(
        suction_pressure_bar,
        discharge_pressure_bar,
        motive_mass,
        suction_mass,
    );
    if mixing_pressure_bar > suction_pressure_bar {
        mixing_pressure_bar = suction_pressure_bar;
    }
    let mixing_pressure = pascals(mixing_pressure_bar * BAR);

    // The motive nozzle: isentropic to the mixing pressure, then to the enthalpy its own
    // efficiency leaves.
    let motive_entropy = motive.entropy()?;
    let motive_mass_enthalpy = motive.h.value / motive.molar_mass()?.value;
    let (_, motive_isentropic_h) = isentropic(motive, mixing_pressure, motive_entropy)?;
    let motive_actual_h = motive_mass_enthalpy
        - setup.motive_nozzle_efficiency * (motive_mass_enthalpy - motive_isentropic_h);
    // The motive's own state at the mixing pressure is the class's *area* calculation -
    // `rhoNozzle` and the Mach number - and neither reaches the discharge state.
    let nozzle_velocity = (2.0 * (motive_mass_enthalpy - motive_actual_h).max(0.0)).sqrt();

    // The suction nozzle, the same route with its own efficiency.
    let suction_entropy = suction.entropy()?;
    let suction_mass_enthalpy = suction.h.value / suction.molar_mass()?.value;
    let (_, suction_isentropic_h) = isentropic(suction, mixing_pressure, suction_entropy)?;
    let suction_actual_h = suction_mass_enthalpy
        - setup.suction_nozzle_efficiency * (suction_mass_enthalpy - suction_isentropic_h);
    let drawn = flash_at_enthalpy(suction, mixing_pressure, suction_actual_h)?;
    let suction_from_enthalpy = (2.0 * (suction_mass_enthalpy - suction_actual_h).max(0.0)).sqrt();
    let drawn_density = drawn.corrected_density()?.max(1.0e-9);
    let design_suction_velocity = suction_from_enthalpy.max(estimate_suction_velocity(
        suction_pressure_bar,
        discharge_pressure_bar,
        drawn_density,
        suction_mass,
    ));

    // The momentum balance, with the mixing chamber's efficiency on the ideal velocity.
    let total_motive_h = motive_actual_h + 0.5 * nozzle_velocity * nozzle_velocity;
    let total_suction_h =
        suction_actual_h + 0.5 * design_suction_velocity * design_suction_velocity;
    let ideal_mixing_velocity =
        (motive_mass * nozzle_velocity + suction_mass * design_suction_velocity) / total_mass;
    let mixing_velocity = setup.mixing_efficiency * ideal_mixing_velocity;
    let mixed_total_h =
        (motive_mass * total_motive_h + suction_mass * total_suction_h) / total_mass;
    let mixed_static_h = mixed_total_h - 0.5 * mixing_velocity * mixing_velocity;

    // The joined fluid at the mixing pressure and that static enthalpy.
    let joined = join(motive, suction)?;
    let _mixing = flash_at_enthalpy(&joined, mixing_pressure, mixed_static_h)?;

    // The diffuser: the recovered energy goes back in as enthalpy, and the pressure rise is
    // taken at the discharge pressure.
    let recovered = setup.diffuser_efficiency * 0.5 * mixing_velocity * mixing_velocity;
    let before_diffuser = mixed_static_h + recovered;
    let at_discharge = flash_at_enthalpy(&joined, setup.discharge_pressure, before_diffuser)?;
    let diffuser_density = at_discharge.corrected_density()?.max(1.0e-9);
    let design_diffuser_velocity = estimate_diffuser_velocity(
        mixing_pressure_bar,
        discharge_pressure_bar,
        diffuser_density,
        total_mass,
    );
    let final_h = before_diffuser - 0.5 * design_diffuser_velocity * design_diffuser_velocity;

    // **Everything the route reached, kept rather than dropped.** The mixing pressure is the
    // class's own estimate and the four velocities are what the efficiencies made of it, so the
    // discharge state is the last line of an argument that was otherwise thrown away.
    Ok(EjectorOutcome {
        outlet: flash_at_enthalpy(&joined, setup.discharge_pressure, final_h)?,
        numbers: Some(EjectorNumbers {
            mixing_pressure,
            motive_nozzle_velocity: meters_per_second(nozzle_velocity),
            suction_nozzle_velocity: meters_per_second(design_suction_velocity),
            mixing_velocity: meters_per_second(mixing_velocity),
            diffuser_velocity: meters_per_second(design_diffuser_velocity),
        }),
    })
}

/// The isentropic flash of one stream to a pressure, as `(temperature, mass enthalpy)`.
fn isentropic(stream: &Stream, p_out: Pressure, s_in: f64) -> Result<(f64, f64)> {
    let (mixture, ideal_gas) = stream.mixture()?;
    let flash = ps_flash(
        &mixture,
        &ideal_gas,
        p_out,
        joules_per_mole_kelvin(s_in),
        &stream.z,
    )?;
    let (h, _) = ph_flash::enthalpy_at(&mixture, &ideal_gas, flash.temperature, p_out, &stream.z)?;
    Ok((flash.temperature.value, h / stream.molar_mass()?.value))
}

/// A stream at a stated **mass** enthalpy and pressure, which is how the class flashes.
///
/// `PHflash(h, "J/kg")` is per kilogram, and this library's flash is per mole, so the mixing
/// and the diffuser steps cross through the mixture's molar mass. The state is a function of
/// `(T, P, z)` either way, and the record this returns is the state's.
fn flash_at_enthalpy(stream: &Stream, p_out: Pressure, h_mass: f64) -> Result<Stream> {
    let (mixture, ideal_gas) = stream.mixture()?;
    let molar_mass = stream.molar_mass()?.value;
    let flash = ph_flash::ph_flash(
        &mixture,
        &ideal_gas,
        p_out,
        joules_per_mole(h_mass * molar_mass),
        &stream.z,
    )?;
    let stage = Stream::from_pt(
        stream.components.clone(),
        stream.z.clone(),
        stream.n,
        p_out,
        flash.temperature,
    )?;
    Ok(stage)
}

/// The two streams as one, by moles - `SystemThermo.addFluid`'s arithmetic.
fn join(motive: &Stream, suction: &Stream) -> Result<Stream> {
    crate::kernels::mixer::mixer(&[motive.clone(), suction.clone()], Some(motive.p))
}

/// The mixing pressure `run` estimates when `setMixingPressure` was never called.
///
/// **In bara, which is the unit the class wrote it in.** The floor is `0.01` bar, not 1000 Pa,
/// and the margin is `3 %` of the suction pressure.
fn estimate_mixing_pressure(
    suction_pressure: f64,
    discharge_pressure: f64,
    motive_mass: f64,
    suction_mass: f64,
) -> f64 {
    if suction_pressure <= 0.0 {
        return discharge_pressure.max(0.0);
    }
    let entrainment_ratio = if motive_mass > 1.0e-9 {
        suction_mass.max(0.0) / motive_mass
    } else {
        1.0
    };
    let clamped_ratio = clamp(entrainment_ratio, 0.0, 5.0);
    let pressure_lift = (discharge_pressure - suction_pressure).max(0.0);
    let pressure_drop = pressure_lift * (0.1 + 0.03 * clamped_ratio);
    let suction_margin = suction_pressure * 0.03;
    let estimated = suction_pressure - pressure_drop.max(suction_margin);
    clamp(
        estimated,
        (0.01f64).max(suction_pressure * 0.4),
        suction_pressure,
    )
}

/// The suction velocity `run` designs to, in m/s - the class's own blend, clamped.
fn estimate_suction_velocity(
    suction_pressure: f64,
    discharge_pressure: f64,
    suction_density: f64,
    suction_mass: f64,
) -> f64 {
    let density = suction_density.max(1.0e-6);
    let available_lift = ((discharge_pressure - suction_pressure) * BAR).max(0.0);
    let target_dynamic = (available_lift * 0.02).max(500.0);
    let baseline = (2.0 * target_dynamic / density).sqrt();
    let volumetric_flow = if suction_mass > 0.0 {
        suction_mass / density
    } else {
        0.0
    };
    let flow_scaling = if volumetric_flow > 0.0 {
        50.0 + 30.0 * (1.0 + volumetric_flow * 5.0).log10()
    } else {
        50.0
    };
    clamp(0.6 * baseline + 0.4 * flow_scaling, 25.0, 120.0)
}

/// The diffuser outlet velocity `run` designs to, in m/s.
fn estimate_diffuser_velocity(
    mixing_pressure: f64,
    discharge_pressure: f64,
    diffuser_density: f64,
    total_mass: f64,
) -> f64 {
    let density = diffuser_density.max(1.0e-6);
    let pressure_recovery = ((discharge_pressure - mixing_pressure) * BAR).max(0.0);
    let target_dynamic = (pressure_recovery * 0.01).max(250.0);
    let baseline = (2.0 * target_dynamic / density).sqrt();
    let volumetric_flow = if total_mass > 0.0 {
        total_mass / density
    } else {
        0.0
    };
    let flow_scaling = if volumetric_flow > 0.0 {
        20.0 + 15.0 * (1.0 + volumetric_flow * 4.0).log10()
    } else {
        25.0
    };
    clamp(0.5 * baseline + 0.5 * flow_scaling, 10.0, 60.0)
}

fn clamp(value: f64, min: f64, max: f64) -> f64 {
    value.max(min).min(max)
}
