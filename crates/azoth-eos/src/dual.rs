//! Second-order forward automatic differentiation in temperature and molar volume.
//!
//! The solid reference equations (argon's Buckingham lattice and para-hydrogen's Vinet
//! curve) are Helmholtz energies `a(T, V)` whose property set is a rearrangement of the
//! derivatives to second order. NeqSim's `thermo.util.solid` computes those with a
//! hand-rolled dual number carrying `a`, `a_T`, `a_V`, `a_TT`, `a_TV` and `a_VV`; this
//! module is that arithmetic, shared rather than re-transcribed per equation.
//!
//! A [`Dual`] is the pair of truncated Taylor coefficients in `T` and `V`. `value` is
//! the scalar, `dt` and `dv` the first derivatives, and `dtt`, `dtv`, `dvv` the second.
//! Every operation below is the chain and product rule written out once, so a
//! transcription error anywhere in a 20-term lattice sum is a test failure here rather
//! than a silently wrong heat capacity.

use std::sync::OnceLock;

/// The scalar value and its derivatives up to second order in `T` and `V`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual {
    /// The value.
    pub value: f64,
    /// `d/dT`.
    pub dt: f64,
    /// `d/dV`.
    pub dv: f64,
    /// `d^2/dT^2`.
    pub dtt: f64,
    /// `d^2/dT dV`.
    pub dtv: f64,
    /// `d^2/dV^2`.
    pub dvv: f64,
}

impl Dual {
    /// A constant independent of `T` and `V`.
    #[must_use]
    pub fn constant(value: f64) -> Self {
        Self {
            value,
            dt: 0.0,
            dv: 0.0,
            dtt: 0.0,
            dtv: 0.0,
            dvv: 0.0,
        }
    }

    /// The temperature variable.
    #[must_use]
    pub fn temperature(value: f64) -> Self {
        Self {
            value,
            dt: 1.0,
            dv: 0.0,
            dtt: 0.0,
            dtv: 0.0,
            dvv: 0.0,
        }
    }

    /// The molar-volume variable.
    #[must_use]
    pub fn volume(value: f64) -> Self {
        Self {
            value,
            dt: 0.0,
            dv: 1.0,
            dtt: 0.0,
            dtv: 0.0,
            dvv: 0.0,
        }
    }

    /// `f(x)` composed through the scalar value, given `f`, `f'` and `f''` at that value.
    #[must_use]
    pub fn apply_unary(self, f: f64, first: f64, second: f64) -> Self {
        Self {
            value: f,
            dt: first * self.dt,
            dv: first * self.dv,
            dtt: second * self.dt * self.dt + first * self.dtt,
            dtv: second * self.dt * self.dv + first * self.dtv,
            dvv: second * self.dv * self.dv + first * self.dvv,
        }
    }

    /// The exponential, `d exp(x)/dx = exp(x)`.
    #[must_use]
    pub fn exp(self) -> Self {
        let result = self.value.exp();
        self.apply_unary(result, result, result)
    }

    /// A constant power, valid while the value is positive.
    #[must_use]
    pub fn powf(self, exponent: f64) -> Self {
        let result = self.value.powf(exponent);
        self.apply_unary(
            result,
            exponent * self.value.powf(exponent - 1.0),
            exponent * (exponent - 1.0) * self.value.powf(exponent - 2.0),
        )
    }

    /// The reciprocal, `1/x`.
    #[must_use]
    pub fn reciprocal(self) -> Self {
        self.apply_unary(
            1.0 / self.value,
            -1.0 / (self.value * self.value),
            2.0 / (self.value * self.value * self.value),
        )
    }

    /// `ln(1 - exp(-x))`, the harmonic-oscillator free energy, composed through the value.
    #[must_use]
    pub fn log_one_minus_exp_negative(self) -> Self {
        let value = self.value;
        if value > 50.0 {
            let exponential = (-value).exp();
            return self.apply_unary((-exponential).ln_1p(), exponential, -exponential);
        }
        let denominator = value.exp_m1();
        self.apply_unary(
            (-(-value).exp_m1()).ln(),
            1.0 / denominator,
            -value.exp() / (denominator * denominator),
        )
    }

    /// The dimensionless Debye Helmholtz function `ln(1 - exp(-x)) - D3(x)` composed
    /// through the reduced temperature `x = theta_D / T`.
    #[must_use]
    pub fn debye_free_energy(self) -> Self {
        let z = self.value;
        let debye = debye_function3(z);
        let (debye_first, debye_second) = if z < 1.0e-3 {
            (
                -1.0 / 8.0 + z / 30.0 - z.powi(3) / 1260.0 + z.powi(5) / 45360.0,
                1.0 / 30.0 - z.powi(2) / 420.0 + z.powi(4) / 9072.0,
            )
        } else {
            let (q, q_first) = if z > 50.0 {
                let exponential = (-z).exp();
                (exponential, -exponential)
            } else {
                let denominator = z.exp_m1();
                let exponential = z.exp();
                (
                    1.0 / denominator,
                    -exponential / (denominator * denominator),
                )
            };
            let first = q - 3.0 * debye / z;
            (first, q_first - 3.0 * first / z + 3.0 * debye / (z * z))
        };

        let (logarithm, logarithm_first, logarithm_second) = if z > 50.0 {
            let exponential = (-z).exp();
            ((-exponential).ln_1p(), exponential, -exponential)
        } else {
            let denominator = z.exp_m1();
            (
                (-(-z).exp_m1()).ln(),
                1.0 / denominator,
                -z.exp() / (denominator * denominator),
            )
        };

        self.apply_unary(
            logarithm - debye,
            logarithm_first - debye_first,
            logarithm_second - debye_second,
        )
    }
}

impl std::ops::Add<Dual> for Dual {
    type Output = Dual;
    fn add(self, other: Dual) -> Dual {
        Dual {
            value: self.value + other.value,
            dt: self.dt + other.dt,
            dv: self.dv + other.dv,
            dtt: self.dtt + other.dtt,
            dtv: self.dtv + other.dtv,
            dvv: self.dvv + other.dvv,
        }
    }
}

impl std::ops::Add<f64> for Dual {
    type Output = Dual;
    fn add(self, scalar: f64) -> Dual {
        Dual {
            value: self.value + scalar,
            ..self
        }
    }
}

impl std::ops::Add<Dual> for f64 {
    type Output = Dual;
    fn add(self, other: Dual) -> Dual {
        other + self
    }
}

impl std::ops::Sub<Dual> for Dual {
    type Output = Dual;
    fn sub(self, other: Dual) -> Dual {
        Dual {
            value: self.value - other.value,
            dt: self.dt - other.dt,
            dv: self.dv - other.dv,
            dtt: self.dtt - other.dtt,
            dtv: self.dtv - other.dtv,
            dvv: self.dvv - other.dvv,
        }
    }
}

impl std::ops::Sub<f64> for Dual {
    type Output = Dual;
    fn sub(self, scalar: f64) -> Dual {
        Dual {
            value: self.value - scalar,
            ..self
        }
    }
}

impl std::ops::Sub<Dual> for f64 {
    type Output = Dual;
    fn sub(self, other: Dual) -> Dual {
        Dual::constant(self) - other
    }
}

impl std::ops::Mul<Dual> for Dual {
    type Output = Dual;
    fn mul(self, other: Dual) -> Dual {
        Dual {
            value: self.value * other.value,
            dt: self.dt * other.value + self.value * other.dt,
            dv: self.dv * other.value + self.value * other.dv,
            dtt: self.dtt * other.value + 2.0 * self.dt * other.dt + self.value * other.dtt,
            dtv: self.dtv * other.value
                + self.dt * other.dv
                + self.dv * other.dt
                + self.value * other.dtv,
            dvv: self.dvv * other.value + 2.0 * self.dv * other.dv + self.value * other.dvv,
        }
    }
}

impl std::ops::Mul<f64> for Dual {
    type Output = Dual;
    fn mul(self, scalar: f64) -> Dual {
        Dual {
            value: self.value * scalar,
            dt: self.dt * scalar,
            dv: self.dv * scalar,
            dtt: self.dtt * scalar,
            dtv: self.dtv * scalar,
            dvv: self.dvv * scalar,
        }
    }
}

impl std::ops::Mul<Dual> for f64 {
    type Output = Dual;
    fn mul(self, other: Dual) -> Dual {
        other * self
    }
}

// Division of a dual by a dual is multiplication by the reciprocal, the standard
// quotient for a truncated power series.
#[allow(clippy::suspicious_arithmetic_impl)]
impl std::ops::Div<Dual> for Dual {
    type Output = Dual;
    fn div(self, other: Dual) -> Dual {
        self * other.reciprocal()
    }
}

impl std::ops::Div<f64> for Dual {
    type Output = Dual;
    fn div(self, scalar: f64) -> Dual {
        Dual {
            value: self.value / scalar,
            dt: self.dt / scalar,
            dv: self.dv / scalar,
            dtt: self.dtt / scalar,
            dtv: self.dtv / scalar,
            dvv: self.dvv / scalar,
        }
    }
}

/// The normalized third-order Debye function, `D3(z) = (1/z^3) * integral_0^z x^3/(e^x - 1) dx`.
///
/// The normalization is the one NeqSim's solid equations use: `D3(0) = 1/3` rather than
/// one. The mid-range value is integrated by 64-point Gauss-Legendre quadrature; the
/// small-`z` series and large-`z` asymptotic are closed forms.
#[must_use]
pub fn debye_function3(z: f64) -> f64 {
    if z < 1.0e-3 {
        return 1.0 / 3.0 - z / 8.0 + z * z / 60.0 - z.powi(4) / 5040.0 + z.powi(6) / 272160.0;
    }
    if z > 50.0 {
        return std::f64::consts::PI.powi(4) / (15.0 * z.powi(3));
    }
    let mut integral = 0.0;
    for (node, weight) in gauss_legendre() {
        let x = 0.5 * z * (node + 1.0);
        integral += weight * x * x * x / x.exp_m1();
    }
    0.5 * z * integral / z.powi(3)
}

/// The molar gas constant NeqSim's phase layer uses for the fugacity, J/(mol.K).
pub const R_FUGACITY: f64 = 8.314_462_1;

/// A solid thermodynamic state at one `(T, P)`, in molar SI units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolidState {
    /// Molar volume, m³/mol.
    pub v: f64,
    /// Molar Helmholtz energy, J/mol.
    pub a: f64,
    /// Molar internal energy, J/mol.
    pub u: f64,
    /// Molar entropy, J/(mol.K).
    pub s: f64,
    /// Molar enthalpy, J/mol.
    pub h: f64,
    /// Molar Gibbs energy, J/mol.
    pub g: f64,
    /// Isochoric heat capacity, J/(mol.K).
    pub cv: f64,
    /// Isobaric heat capacity, J/(mol.K).
    pub cp: f64,
    /// Natural logarithm of the fugacity coefficient, on the 1 bara standard state.
    pub ln_phi: f64,
}

/// Assemble the solid state from a Helmholtz derivative at a temperature, pressure and
/// molar volume. The volume root is the caller's: the pressure here is the one implied by
/// `-a_V` at that volume.
#[must_use]
pub fn state_from_helmholtz(
    temperature: f64,
    pressure_pa: f64,
    molar_volume: f64,
    helmholtz: Dual,
) -> SolidState {
    let entropy = -helmholtz.dt;
    let internal_energy = helmholtz.value + temperature * entropy;
    let gibbs = helmholtz.value + pressure_pa * molar_volume;
    let cv = -temperature * helmholtz.dtt;
    SolidState {
        v: molar_volume,
        a: helmholtz.value,
        u: internal_energy,
        s: entropy,
        h: internal_energy + pressure_pa * molar_volume,
        g: gibbs,
        cv,
        cp: cv + temperature * helmholtz.dtv * helmholtz.dtv / helmholtz.dvv,
        ln_phi: gibbs / (R_FUGACITY * temperature) - (pressure_pa / 1.0e5).ln(),
    }
}

const QUADRATURE_ORDER: usize = 64;

/// The 64-point Gauss-Legendre nodes and weights, computed once.
fn gauss_legendre() -> &'static [(f64, f64)] {
    static NODES: OnceLock<Vec<(f64, f64)>> = OnceLock::new();
    NODES.get_or_init(|| {
        let mut nodes = vec![(0.0, 0.0); QUADRATURE_ORDER];
        let midpoint = QUADRATURE_ORDER.div_ceil(2);
        for i in 0..midpoint {
            let mut root =
                ((i as f64 + 0.75) * std::f64::consts::PI / (QUADRATURE_ORDER as f64 + 0.5)).cos();
            let derivative = loop {
                let previous = root;
                let mut p0 = 1.0;
                let mut p1 = root;
                for order in 2..=QUADRATURE_ORDER {
                    let polynomial = ((2.0 * order as f64 - 1.0) * root * p1
                        - (order as f64 - 1.0) * p0)
                        / order as f64;
                    p0 = p1;
                    p1 = polynomial;
                }
                let derivative = QUADRATURE_ORDER as f64 * (root * p1 - p0) / (root * root - 1.0);
                root = previous - p1 / derivative;
                if (root - previous).abs() <= 1.0e-15 {
                    break derivative;
                }
            };
            let weight = 2.0 / ((1.0 - root * root) * derivative * derivative);
            nodes[i] = (-root, weight);
            nodes[QUADRATURE_ORDER - 1 - i] = (root, weight);
        }
        nodes
    })
}
