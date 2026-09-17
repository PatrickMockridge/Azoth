//! The GERG-2008 multi-fluid Helmholtz equation of state, as NeqSim's `GERG2008`
//! carries it.
//!
//! GERG-2008 (Kunz and Wagner 2012) is the 21-component wide-range mixture model for
//! natural gas: a mixture Helmholtz energy in reduced variables `tau = Tr/T`,
//! `delta = D/Dr` with the reducing temperature and density `(Tr, Dr)` a composition
//! average over binary interaction parameters, plus a pairwise departure contribution.
//! The coefficients live in [`crate::gerg2008_data`]; this module is the arithmetic.
//!
//! Native units are NeqSim's: temperature in K, density in mol/L, pressure in kPa, and
//! `R = 8.314472` J/(mol.K). The model layer converts to the molar SI convention.

// The 1-indexed loops mirror NeqSim's Fortran-style array indices so the port reads
// against the source; an iterator rewrite would obscure that correspondence.
#![allow(clippy::needless_range_loop)]

use crate::gerg2008_data as data;

/// The gas constant of the equation, J/(mol.K).
pub const R: f64 = data::R;

/// The number of components.
pub const NCOMP: usize = data::NCOMP;

/// The component names, indexed as NeqSim's GERG-2008 does (index 0 is unused).
pub const COMPONENT_NAMES: [&str; NCOMP + 1] = [
    "",
    "methane",
    "nitrogen",
    "CO2",
    "ethane",
    "propane",
    "i-butane",
    "n-butane",
    "i-pentane",
    "n-pentane",
    "n-hexane",
    "n-heptane",
    "n-octane",
    "n-nonane",
    "nC10",
    "hydrogen",
    "oxygen",
    "CO",
    "water",
    "H2S",
    "helium",
    "argon",
];

const EPSILON: f64 = 1e-15;

/// The component index for a name, or `None` if the name is not a GERG-2008 component.
#[must_use]
pub fn component_index(name: &str) -> Option<usize> {
    COMPONENT_NAMES.iter().position(|c| *c == name)
}

/// The molar thermodynamic state of one GERG-2008 mixture at `(T, D, x)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Properties {
    /// Pressure, kPa.
    pub pressure_kpa: f64,
    /// Compressibility factor.
    pub z: f64,
    /// Internal energy, J/mol.
    pub u: f64,
    /// Enthalpy, J/mol.
    pub h: f64,
    /// Entropy, J/(mol.K).
    pub s: f64,
    /// Isochoric heat capacity, J/(mol.K).
    pub cv: f64,
    /// Isobaric heat capacity, J/(mol.K).
    pub cp: f64,
    /// Gibbs energy, J/mol.
    pub g: f64,
    /// Speed of sound, m/s.
    pub w: f64,
}

/// The reducing temperature and density `(Tr, Dr)` of a composition.
///
/// `x` is 1-indexed, length `NCOMP + 1`, summing to one.
#[must_use]
pub fn reducing_parameters(x: &[f64]) -> (f64, f64) {
    let mut tr = 0.0;
    let mut vr = 0.0;
    for i in 1..=NCOMP {
        if x[i] > EPSILON {
            let mut f = 1.0;
            for j in i..=NCOMP {
                if x[j] > EPSILON {
                    let xij = f * (x[i] * x[j]) * (x[i] + x[j]);
                    vr += xij * data::GVIJ[i][j] / (data::BVIJ[i][j] * x[i] + x[j]);
                    tr += xij * data::GTIJ[i][j] / (data::BTIJ[i][j] * x[i] + x[j]);
                    f = 2.0;
                }
            }
        }
    }
    let dr = if vr > EPSILON { 1.0 / vr } else { 0.0 };
    (tr, dr)
}

/// The ideal-gas Helmholtz energy and its first two `tau` derivatives, as `a0 = [alpha0,
/// tau*d(alpha0)/dtau, tau^2*d^2(alpha0)/dtau^2]`.
#[must_use]
pub fn alpha0(t: f64, d: f64, x: &[f64]) -> [f64; 3] {
    let log_d = if d > EPSILON { d.ln() } else { EPSILON.ln() };
    let log_t = t.ln();
    let mut a0 = [0.0; 3];
    for i in 1..=NCOMP {
        if x[i] > EPSILON {
            let log_xd = log_d + x[i].ln();
            let mut sum_hyp0 = 0.0;
            let mut sum_hyp1 = 0.0;
            let mut sum_hyp2 = 0.0;
            for j in 4..=7 {
                if data::TH0I[i][j] > EPSILON {
                    let th0t = data::TH0I[i][j] / t;
                    let ep = th0t.exp();
                    let em = 1.0 / ep;
                    let hsn = (ep - em) / 2.0;
                    let hcn = (ep + em) / 2.0;
                    if j == 4 || j == 6 {
                        let log_hyp = hsn.abs().ln();
                        sum_hyp0 += data::N0I[i][j] * log_hyp;
                        sum_hyp1 += data::N0I[i][j] * th0t * hcn / hsn;
                        sum_hyp2 += data::N0I[i][j] * (th0t / hsn) * (th0t / hsn);
                    } else {
                        let log_hyp = hcn.abs().ln();
                        sum_hyp0 -= data::N0I[i][j] * log_hyp;
                        sum_hyp1 -= data::N0I[i][j] * th0t * hsn / hcn;
                        sum_hyp2 += data::N0I[i][j] * (th0t / hcn) * (th0t / hcn);
                    }
                }
            }
            a0[0] += x[i]
                * (log_xd + data::N0I[i][1] + data::N0I[i][2] / t - data::N0I[i][3] * log_t
                    + sum_hyp0);
            a0[1] += x[i] * (data::N0I[i][3] + data::N0I[i][2] / t + sum_hyp1);
            a0[2] += -x[i] * (data::N0I[i][3] + sum_hyp2);
        }
    }
    a0
}

/// The tau-dependent coefficient parts, `taup[i][k] = noik * tau^toik`, with the GERG-2008
/// short-form sharing of the propane exponents.
fn t_terms(lntau: f64, x: &[f64], taup: &mut [[f64; 25]; NCOMP + 1]) {
    let mut taup0 = [0.0; 13];
    for k in 1..=(data::KPOL[5] + data::KEXP[5]) as usize {
        taup0[k] = (data::TOIK[5][k] * lntau).exp();
    }
    for i in 1..=NCOMP {
        if x[i] > EPSILON {
            let short = i > 4 && i != 15 && i != 18 && i != 20;
            for k in 1..=(data::KPOL[i] + data::KEXP[i]) as usize {
                if short {
                    taup[i][k] = data::NOIK[i][k] * taup0[k];
                } else {
                    taup[i][k] = data::NOIK[i][k] * (data::TOIK[i][k] * lntau).exp();
                }
            }
        }
    }
}

/// The residual Helmholtz derivatives, indexed `ar[tau-order][delta-order]`.
///
/// `itau` computes the tau derivatives (`ar[1][..]`, `ar[2][0]`) and the third delta
/// derivative; the pressure-only call sets it false.
#[must_use]
pub fn alphar(itau: bool, t: f64, d: f64, x: &[f64]) -> [[f64; 4]; 4] {
    let mut ar = [[0.0; 4]; 4];
    let (tr, dr) = reducing_parameters(x);
    let del = d / dr;
    let tau = tr / t;
    let lntau = tau.ln();

    let mut delp = [0.0; 8];
    let mut expd = [0.0; 8];
    delp[1] = del;
    expd[1] = (-delp[1]).exp();
    for i in 2..=7 {
        delp[i] = delp[i - 1] * del;
        expd[i] = (-delp[i]).exp();
    }

    let mut taup = [[0.0; 25]; NCOMP + 1];
    t_terms(lntau, x, &mut taup);

    // Pure-fluid contributions.
    for i in 1..=NCOMP {
        if x[i] <= EPSILON {
            continue;
        }
        for k in 1..=data::KPOL[i] as usize {
            let ndt = x[i] * delp[data::DOIK[i][k] as usize] * taup[i][k];
            let ndtd = ndt * data::DOIK[i][k] as f64;
            ar[0][1] += ndtd;
            ar[0][2] += ndtd * (data::DOIK[i][k] - 1) as f64;
            if itau {
                let ndtt = ndt * data::TOIK[i][k];
                ar[0][0] += ndt;
                ar[1][0] += ndtt;
                ar[2][0] += ndtt * (data::TOIK[i][k] - 1.0);
                ar[1][1] += ndtt * data::DOIK[i][k] as f64;
                ar[1][2] += ndtt * data::DOIK[i][k] as f64 * (data::DOIK[i][k] - 1) as f64;
                ar[0][3] += ndtd * (data::DOIK[i][k] - 1) as f64 * (data::DOIK[i][k] - 2) as f64;
            }
        }
        for k in (1 + data::KPOL[i]) as usize..=(data::KPOL[i] + data::KEXP[i]) as usize {
            let ndt = x[i]
                * delp[data::DOIK[i][k] as usize]
                * taup[i][k]
                * expd[data::COIK[i][k] as usize];
            let ex = data::COIK[i][k] as f64 * delp[data::COIK[i][k] as usize];
            let ex2 = data::DOIK[i][k] as f64 - ex;
            let ex3 = ex2 * (ex2 - 1.0);
            ar[0][1] += ndt * ex2;
            ar[0][2] += ndt * (ex3 - data::COIK[i][k] as f64 * ex);
            if itau {
                let ndtt = ndt * data::TOIK[i][k];
                ar[0][0] += ndt;
                ar[1][0] += ndtt;
                ar[2][0] += ndtt * (data::TOIK[i][k] - 1.0);
                ar[1][1] += ndtt * ex2;
                ar[1][2] += ndtt * (ex3 - data::COIK[i][k] as f64 * ex);
                ar[0][3] += ndt
                    * (ex3 * (ex2 - 2.0)
                        - ex * (3.0 * ex2 - 3.0 + data::COIK[i][k] as f64)
                            * data::COIK[i][k] as f64);
            }
        }
    }

    // Mixture departure contributions.
    let mut taupijk = [[0.0; 13]; 15];
    for i in 1..NCOMP {
        if x[i] <= EPSILON {
            continue;
        }
        for j in (i + 1)..=NCOMP {
            if x[j] <= EPSILON {
                continue;
            }
            let mn = data::MNUMB[i][j];
            if mn < 0 {
                continue;
            }
            let mn = mn as usize;
            let xijf = x[i] * x[j] * data::FIJ[i][j];
            for k in 1..=data::KPOLIJ[mn] as usize {
                taupijk[mn][k] = data::NIJK[mn][k] * (data::TIJK[mn][k] * lntau).exp();
                let ndt = xijf * delp[data::DIJK[mn][k] as usize] * taupijk[mn][k];
                let ndtd = ndt * data::DIJK[mn][k] as f64;
                ar[0][1] += ndtd;
                ar[0][2] += ndtd * (data::DIJK[mn][k] - 1) as f64;
                if itau {
                    let ndtt = ndt * data::TIJK[mn][k];
                    ar[0][0] += ndt;
                    ar[1][0] += ndtt;
                    ar[2][0] += ndtt * (data::TIJK[mn][k] - 1.0);
                    ar[1][1] += ndtt * data::DIJK[mn][k] as f64;
                    ar[1][2] += ndtt * data::DIJK[mn][k] as f64 * (data::DIJK[mn][k] - 1) as f64;
                    ar[0][3] +=
                        ndtd * (data::DIJK[mn][k] - 1) as f64 * (data::DIJK[mn][k] - 2) as f64;
                }
            }
            for k in
                (1 + data::KPOLIJ[mn]) as usize..=(data::KPOLIJ[mn] + data::KEXPIJ[mn]) as usize
            {
                let cij0 = data::CIJK[mn][k] * delp[2];
                let eij0 = data::EIJK[mn][k] * del;
                let ndt = xijf
                    * data::NIJK[mn][k]
                    * delp[data::DIJK[mn][k] as usize]
                    * (cij0 + eij0 + data::GIJK[mn][k] + data::TIJK[mn][k] * lntau).exp();
                let ex = data::DIJK[mn][k] as f64 + 2.0 * cij0 + eij0;
                let ex2 = ex * ex - data::DIJK[mn][k] as f64 + 2.0 * cij0;
                ar[0][1] += ndt * ex;
                ar[0][2] += ndt * ex2;
                if itau {
                    let ndtt = ndt * data::TIJK[mn][k];
                    ar[0][0] += ndt;
                    ar[1][0] += ndtt;
                    ar[2][0] += ndtt * (data::TIJK[mn][k] - 1.0);
                    ar[1][1] += ndtt * ex;
                    ar[1][2] += ndtt * ex2;
                    ar[0][3] += ndt
                        * (ex * (ex2 - 2.0 * (data::DIJK[mn][k] as f64 - 2.0 * cij0))
                            + 2.0 * data::DIJK[mn][k] as f64);
                }
            }
        }
    }

    ar
}

/// The pressure (kPa) and its density derivative `dPdDsave` at `(T, D, x)`.
fn pressure(t: f64, d: f64, x: &[f64]) -> (f64, f64) {
    let ar = alphar(false, t, d, x);
    let z = 1.0 + ar[0][1];
    let p = d * R * t * z;
    let dpdd = R * t * (1.0 + 2.0 * ar[0][1] + ar[0][2]);
    (p, dpdd)
}

/// The molar mass of a composition, g/mol.
#[must_use]
pub fn molar_mass(x: &[f64]) -> f64 {
    let mut mm = 0.0;
    for i in 1..=NCOMP {
        mm += x[i] * data::MM[i];
    }
    mm
}

/// The full property set at a temperature (K), density (mol/L) and composition.
#[must_use]
pub fn properties(t: f64, d: f64, x: &[f64]) -> Properties {
    let mm = molar_mass(x);
    let a0 = alpha0(t, d, x);
    let ar = alphar(true, t, d, x);

    let rt = R * t;
    let z = 1.0 + ar[0][1];
    let pressure_kpa = d * rt * z;
    let dpdd = rt * (1.0 + 2.0 * ar[0][1] + ar[0][2]);
    let dpdt = d * R * (1.0 + ar[0][1] - ar[1][1]);

    let u = rt * (a0[1] + ar[1][0]);
    let h = rt * (1.0 + ar[0][1] + a0[1] + ar[1][0]);
    let s = R * (a0[1] + ar[1][0] - a0[0] - ar[0][0]);
    let cv = -R * (a0[2] + ar[2][0]);
    let g = rt * (1.0 + ar[0][1] + a0[0] + ar[0][0]);

    let cp = if d > EPSILON {
        cv + t * (dpdt / d) * (dpdt / d) / dpdd
    } else {
        cv + R
    };

    let w = if cv > EPSILON {
        (1000.0 * cp / cv * dpdd / mm).max(0.0).sqrt()
    } else {
        0.0
    };

    Properties {
        pressure_kpa,
        z,
        u,
        h,
        s,
        cv,
        cp,
        g,
        w,
    }
}

/// Solve `P(D) = pressure_kpa` for the molar density (mol/L) by Newton iteration in
/// `log(1/D)`, with the GERG-2008 restart strategy.
#[must_use]
pub fn solve_density(t: f64, pressure_kpa: f64, x: &[f64]) -> f64 {
    let (tcx, dcx) = pseudo_critical_point(x);
    let mut d = pressure_kpa / R / t;
    let plog = pressure_kpa.ln();
    let mut vlog = -d.ln();
    let tolr = 1e-7;
    let mut n_fail = 0;
    let mut i_fail = 0;

    for it in 1..=50 {
        if !(-7.0..=100.0).contains(&vlog) || matches!(it, 20 | 30 | 40) || i_fail == 1 {
            i_fail = 0;
            if n_fail > 2 {
                return pressure_kpa / R / t;
            }
            n_fail += 1;
            d = match n_fail {
                1 => dcx * 3.0,
                2 => dcx * 2.5,
                _ => dcx * 2.0,
            };
            vlog = -d.ln();
        }
        d = (-vlog).exp();
        let (p2, dpdd) = pressure(t, d, x);
        if dpdd < EPSILON || p2 < EPSILON {
            let mut vinc = if d > dcx { -0.1 } else { 0.1 };
            if it > 5 {
                vinc /= 2.0;
            }
            if it > 10 && it < 20 {
                vinc /= 5.0;
            }
            vlog += vinc;
        } else {
            let dpdlv = -d * dpdd;
            let vdiff = (p2.ln() - plog) * p2 / dpdlv;
            vlog += -vdiff;
            if vdiff.abs() < tolr {
                if dpdd < 0.0 {
                    i_fail = 1;
                } else {
                    let _ = tcx;
                    return (-vlog).exp();
                }
            }
        }
    }
    pressure_kpa / R / t
}

/// The pseudo-critical temperature and density, mole-fraction averages.
fn pseudo_critical_point(x: &[f64]) -> (f64, f64) {
    let mut tcx = 0.0;
    let mut vcx = 0.0;
    for i in 1..=NCOMP {
        tcx += x[i] * data::TC[i];
        vcx += x[i] / data::DC[i];
    }
    let dcx = if vcx > EPSILON { 1.0 / vcx } else { 0.0 };
    (tcx, dcx)
}
