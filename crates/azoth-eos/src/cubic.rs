//! The cubic equation of state's shape: the constants that fix one cubic.
//!
//! The general cubic
//!
//! ```text
//! P = R T/(v - b) - a/((v + delta1 b)(v + delta2 b))
//! ```
//!
//! is fixed by four numbers - `omega_a`, `omega_b`, `delta1`, `delta2` - and the three
//! models this library ports are three settings of them. The temperature dependence of
//! `a` is a separate axis, in [`crate::alpha_term`], because NeqSim decouples the two:
//! a component owns the shape here, and a per-component term owns `alpha(T)`.

/// A cubic equation of state's shape.
///
/// The whole model layer reads its geometry through one of these rather than inlining
/// `sqrt(2)` the way the Peng-Robinson-only code once did, so a second cubic is a new
/// variant rather than a new branch in every fugacity and Helmholtz expression.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Cubic {
    /// Peng-Robinson: `omega = (0.45724333333, 0.077803333)`, `delta = (1 + sqrt(2), 1 - sqrt(2))`.
    ///
    /// The `omega` pair is NeqSim's, not the paper's - see
    /// [`crate::pr_alpha_ab::OMEGA_A`] for the argument.
    #[default]
    Pr,
    /// Soave-Redlich-Kwong: `omega = (1/(9(2^(1/3) - 1)), (2^(1/3) - 1)/3)`, `delta = (1, 0)`.
    ///
    /// NeqSim's `ComponentSrk` computes the pair from `Math.pow(2.0, 1.0/3.0)` at
    /// construction rather than carrying literals; the two decimals below are that
    /// expression to full double precision, and the test
    /// `a_soave_redlich_kwong_recomputes_its_omegas` re-derives them.
    Srk,
    /// Redlich-Kwong: the original, with the same `omega` and `delta` as Soave's but
    /// `alpha = 1/sqrt(Tr)` instead of Soave's correlation. The shape is identical to
    /// [`Cubic::Srk`]; only the alpha term in [`crate::alpha_term`] differs.
    Rk,
    /// Twu-Sim-Tassone: Soave's attraction and repulsion constants with Peng-Robinson's
    /// geometry, and Twu's alpha. NeqSim's `ComponentTST`.
    ///
    /// `omega = (0.427481, 0.086641)` - Soave's constants rounded to six decimals, as
    /// NeqSim writes them - and `delta = (1 + sqrt(2), 1 - sqrt(2))`. The root finder
    /// is shared with Peng-Robinson; only the constants here and the Twu alpha differ.
    Tst,
}

impl Cubic {
    /// `omega_a`, the attraction constant.
    #[must_use]
    pub const fn omega_a(self) -> f64 {
        match self {
            Cubic::Pr => 0.45724333333,
            Cubic::Srk => 0.4274802335403413,
            Cubic::Rk => 0.4274802335403413,
            Cubic::Tst => 0.427481,
        }
    }

    /// `omega_b`, the repulsion constant.
    #[must_use]
    pub const fn omega_b(self) -> f64 {
        match self {
            Cubic::Pr => 0.077803333,
            Cubic::Srk => 0.08664034996495773,
            Cubic::Rk => 0.08664034996495773,
            Cubic::Tst => 0.086641,
        }
    }

    /// `delta1`, the larger geometry delta.
    #[must_use]
    pub const fn delta1(self) -> f64 {
        match self {
            Cubic::Pr => 1.0 + std::f64::consts::SQRT_2,
            Cubic::Srk => 1.0,
            Cubic::Rk => 1.0,
            Cubic::Tst => 1.0 + std::f64::consts::SQRT_2,
        }
    }

    /// `delta2`, the smaller geometry delta.
    #[must_use]
    pub const fn delta2(self) -> f64 {
        match self {
            Cubic::Pr => 1.0 - std::f64::consts::SQRT_2,
            Cubic::Srk => 0.0,
            Cubic::Rk => 0.0,
            Cubic::Tst => 1.0 - std::f64::consts::SQRT_2,
        }
    }

    /// `delta1 - delta2`, the `2 sqrt(2)` of the Peng-Robinson forms.
    #[must_use]
    pub const fn delta_diff(self) -> f64 {
        match self {
            Cubic::Pr => 2.0 * std::f64::consts::SQRT_2,
            Cubic::Srk => 1.0,
            Cubic::Rk => 1.0,
            Cubic::Tst => 2.0 * std::f64::consts::SQRT_2,
        }
    }

    /// `delta1 + delta2`.
    #[must_use]
    pub const fn delta_sum(self) -> f64 {
        match self {
            Cubic::Pr => 2.0,
            Cubic::Srk => 1.0,
            Cubic::Rk => 1.0,
            Cubic::Tst => 2.0,
        }
    }

    /// `delta1 * delta2`.
    #[must_use]
    pub const fn delta_prod(self) -> f64 {
        match self {
            Cubic::Pr => -1.0,
            Cubic::Srk => 0.0,
            Cubic::Rk => 0.0,
            Cubic::Tst => -1.0,
        }
    }

    /// `(delta1 - delta2)/2`, the `sqrt(2)` that the Helmholtz terms carry bare.
    #[must_use]
    pub fn half_delta_diff(self) -> f64 {
        self.delta_diff() / 2.0
    }

    /// `ln((z + delta1 B)/(z + delta2 B))`, the "I" term of the fugacity coefficient.
    #[must_use]
    pub fn i_term(self, z: f64, b: f64) -> f64 {
        ((z + self.delta1() * b) / (z + self.delta2() * b)).ln()
    }

    /// `A/((delta1 - delta2) B)`, the coefficient the fugacity term is scaled by.
    #[must_use]
    pub fn coefficient(self, a: f64, b: f64) -> f64 {
        a / (self.delta_diff() * b)
    }

    /// `z - 1` as the equation of state writes it: `B/(z - B) - A z/((z + d1 B)(z + d2 B))`.
    ///
    /// **Equal to `z - 1` at every root of this cubic, and to nothing else anywhere else.**
    /// The fugacity coefficient's first term is `(B_i/B)` times this, because that term is a
    /// derivative of the Helmholtz energy *at the volume it is given*: a plain cubic sits at
    /// its own root and `z - 1` will do, an associating mixture's volume is moved by the
    /// association's pressure and it will not - by `0.0328` in `ln phi_0` at 356 K and
    /// 1 bara, which is the whole of the divergence that was first read as NeqSim's.
    #[must_use]
    pub fn eos_z_minus_one(self, z: f64, a: f64, b: f64) -> f64 {
        self.eos_z_minus_one_partials(z, a, b).0
    }

    /// [`Self::eos_z_minus_one`], with its `d/dz`, `d/dA` and `d/dB`.
    #[must_use]
    pub fn eos_z_minus_one_partials(self, z: f64, a: f64, b: f64) -> (f64, f64, f64, f64) {
        let (d1, d2) = (self.delta1(), self.delta2());
        let (p1, p2) = (z + d1 * b, z + d2 * b);
        let d = p1 * p2;
        let value = b / (z - b) - a * z / d;
        let dz = -b / ((z - b) * (z - b)) - a * (d - z * (p1 + p2)) / (d * d);
        let da = -z / d;
        let db = z / ((z - b) * (z - b)) + a * z * (d1 * p2 + d2 * p1) / (d * d);
        (value, dz, da, db)
    }

    /// The Huron-Vidal constant `ln((1+delta1)/(1+delta2)) / (delta1 - delta2)`.
    ///
    /// The geometry factor that turns the excess Gibbs energy at infinite pressure into
    /// the attraction parameter. `ln(2)` for Soave, `ln((2+sqrt2)/(2-sqrt2))/(2 sqrt2)`
    /// for Peng-Robinson.
    #[must_use]
    pub fn hv_constant(self) -> f64 {
        ((1.0 + self.delta1()) / (1.0 + self.delta2())).ln() / self.delta_diff()
    }

    /// The monic cubic `z**3 + c2 z**2 + c1 z + c0` in the reduced parameters.
    #[must_use]
    pub fn z_coefficients(self, a: f64, b: f64) -> (f64, f64, f64) {
        let ds = self.delta_sum();
        let dp = self.delta_prod();
        let c2 = (ds - 1.0) * b - 1.0;
        let c1 = dp * b * b - ds * b * b - ds * b + a;
        let c0 = -dp * b * b * b - dp * b * b - a * b;
        (c2, c1, c0)
    }

    /// `dF/dz` of the z-cubic, at a root `z`.
    #[must_use]
    pub fn df_dz(self, z: f64, a: f64, b: f64) -> f64 {
        let (c2, c1, _) = self.z_coefficients(a, b);
        3.0 * z * z + 2.0 * c2 * z + c1
    }

    /// `dF/dA` of the z-cubic, at a root `z`.
    ///
    /// `A` appears in `c1` with coefficient one and in `c0` as `-A B`, so this is
    /// `z - B` for every cubic - which is what makes the implicit `dZ/dA` a single
    /// division below.
    #[must_use]
    pub fn df_da(self, z: f64, b: f64) -> f64 {
        z - b
    }

    /// `dF/dB` of the z-cubic, at a root `z`.
    #[must_use]
    pub fn df_db(self, z: f64, a: f64, b: f64) -> f64 {
        let ds = self.delta_sum();
        let dp = self.delta_prod();
        // From `c2 = (delta1 + delta2 - 1) B - 1`, `c1` and `c0` in `z_coefficients`.
        let dc2_db = ds - 1.0;
        let dc1_db = 2.0 * (dp - ds) * b - ds;
        let dc0_db = -3.0 * dp * b * b - 2.0 * dp * b - a;
        dc2_db * z * z + dc1_db * z + dc0_db
    }

    /// `T * dF/dT` of the z-cubic, with `t_da = T * da/dT` and `t_db = T * db/dT`.
    #[must_use]
    pub fn t_dfdt(self, z: f64, a: f64, b: f64, t_da: f64, t_db: f64) -> f64 {
        let ds = self.delta_sum();
        let dp = self.delta_prod();
        // `T * d(c2)/dT = (delta1 + delta2 - 1) t_db`.
        let t_dc2 = (ds - 1.0) * t_db;
        // `T * d(c1)/dT = (2 B (delta1 delta2 - delta1 - delta2) - (delta1 + delta2)) t_db + t_da`.
        let t_dc1 = (2.0 * b * (dp - ds) - ds) * t_db + t_da;
        // `T * d(c0)/dT = (-3 delta1 delta2 B**2 - 2 delta1 delta2 B - A) t_db - t_da B`.
        let t_dc0 = (-3.0 * dp * b * b - 2.0 * dp * b - a) * t_db - t_da * b;
        t_dc2 * z * z + t_dc1 * z + t_dc0
    }

    /// `G(b) = ln((1 + delta1 b)/(1 + delta2 b))`, the Helmholtz geometry term.
    #[must_use]
    pub fn helmholtz_g(self, b: f64) -> f64 {
        ((1.0 + self.delta1() * b) / (1.0 + self.delta2() * b)).ln()
    }

    /// `q(b) = (1 + delta1 b)(1 + delta2 b)`, the denominator `G'` and `G''` share.
    #[must_use]
    pub fn helmholtz_q(self, b: f64) -> f64 {
        1.0 + self.delta_sum() * b + self.delta_prod() * b * b
    }

    /// `G'(b) = (delta1 - delta2)/q`.
    #[must_use]
    pub fn helmholtz_g_prime(self, b: f64) -> f64 {
        self.delta_diff() / self.helmholtz_q(b)
    }

    /// `G''(b) = -(delta1 - delta2)(delta1 + delta2 + 2 delta1 delta2 b)/q**2`.
    #[must_use]
    pub fn helmholtz_g_second(self, b: f64) -> f64 {
        let q = self.helmholtz_q(b);
        -self.delta_diff() * (self.delta_sum() + 2.0 * self.delta_prod() * b) / (q * q)
    }

    /// The short name that crosses the Python boundary, and the spelling a keycard
    /// would use.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Cubic::Pr => "pr",
            Cubic::Srk => "srk",
            Cubic::Rk => "rk",
            Cubic::Tst => "tst",
        }
    }
}

impl std::str::FromStr for Cubic {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "pr" => Ok(Cubic::Pr),
            "srk" => Ok(Cubic::Srk),
            "rk" => Ok(Cubic::Rk),
            "tst" => Ok(Cubic::Tst),
            other => Err(format!(
                "unknown cubic `{other}`; expected `pr`, `srk`, `rk` or `tst`"
            )),
        }
    }
}
