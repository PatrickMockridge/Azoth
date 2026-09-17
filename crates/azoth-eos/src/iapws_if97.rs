//! The IAPWS-IF97 industrial formulation for water.
//!
//! The Gibbs-energy form of the steam tables: Region 1 (subcooled liquid) and Region 2
//! (superheated vapour) as dimensionless Gibbs functions in the reduced variables
//! `pi = p/p*` and `tau = T*/T`, with Region 4 the saturation curve that selects between
//! them, as NeqSim's `thermo.util.steam.Iapws_if97` carries it. Regions 3 (near-critical)
//! and 5 (high-temperature) are not implemented in NeqSim and are not here.
//!
//! IF97 is mass-based: `R = 0.461526` kJ/(kg.K), pressures in MPa, and the properties
//! specific (kJ/kg, m³/kg). The mixture layer converts to the molar SI convention.

/// The specific gas constant for water, kJ/(kg.K).
pub const R: f64 = 0.461526;

/// The molar mass of water, kg/mol, the single constant that turns a specific property
/// into a molar one.
pub const MOLAR_MASS: f64 = 0.018_015_28;

// Region 1 coefficients (34 terms): `gamma1 = sum N1[i] (7.1 - pi)^I1[i] (tau - 1.222)^J1[i]`.
const I1: [i32; 34] = [
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29,
    30, 31, 32,
];
const J1: [i32; 34] = [
    -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11,
    -6, -29, -31, -38, -39, -40, -41,
];
const N1: [f64; 34] = [
    0.146_329_712_131_67,
    -0.845_481_871_691_14,
    -3.756_360_367_204,
    3.385_516_916_838_5,
    -0.957_919_633_878_72,
    0.157_720_385_132_28,
    -0.016_616_417_199_501,
    0.000_812_146_299_835_68,
    0.000_283_190_801_238_04,
    -0.000_607_063_015_658_74,
    -0.018_990_068_218_419,
    -0.032_529_748_770_505,
    -0.021_841_717_175_414,
    -0.000_052_838_357_969_93,
    -0.000_471_843_210_732_67,
    -0.000_300_017_807_930_26,
    0.000_047_661_393_906_987,
    -0.000_004_414_184_533_084_6,
    -7.269_499_629_759_4e-16,
    -0.000_031_679_644_845_054,
    -0.000_002_827_079_798_531_2,
    -8.520_512_812_010_3e-10,
    -0.000_002_242_528_190_8,
    -6.517_122_289_560_1e-7,
    -1.434_172_993_792_4e-13,
    -4.051_699_686_011_7e-7,
    -1.273_430_174_164_1e-9,
    -1.742_487_123_063_4e-10,
    -6.876_213_129_553_1e-19,
    1.447_830_782_852_1e-20,
    2.633_578_166_279_5e-23,
    -1.194_762_264_007_1e-23,
    1.822_809_458_140_4e-24,
    -9.353_708_729_245_8e-26,
];

// Region 2 ideal part (9 terms): `gamma0 = ln(pi) + sum N0[i] tau^J0[i]`.
const J0: [i32; 9] = [0, 1, -5, -4, -3, -2, -1, 2, 3];
const N0: [f64; 9] = [
    -9.692_768_650_021_7,
    10.086_655_968_018,
    -0.005_608_791_128_302,
    0.071_452_738_081_455,
    -0.407_104_982_239_28,
    1.424_081_917_144_4,
    -4.383_951_131_945,
    -0.284_086_324_607_72,
    0.021_268_463_753_307,
];

// Region 2 residual part (43 terms): `gammar = sum NR[i] pi^IR[i] (tau - 0.5)^JR[i]`.
const IR: [i32; 43] = [
    1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10,
    16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24,
];
const JR: [i32; 43] = [
    0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4,
    10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58,
];
const NR: [f64; 43] = [
    -0.001_773_174_247_321_3,
    -0.017_834_862_292_358,
    -0.045_996_013_696_365,
    -0.057_581_259_083_432,
    -0.050_325_278_727_93,
    -0.000_033_032_641_670_203,
    -0.000_189_489_875_163_15,
    -0.003_939_277_724_335_5,
    -0.043_797_295_650_573,
    -0.000_026_674_547_914_087,
    2.048_173_769_230_9e-8,
    4.387_066_728_443_5e-7,
    -0.000_032_277_677_238_57,
    -0.001_503_392_454_214_8,
    -0.040_668_253_562_649,
    -7.884_730_955_936_7e-10,
    1.279_071_785_228_5e-8,
    4.822_537_271_850_7e-7,
    0.000_002_292_207_633_766_1,
    -1.671_476_645_106_1e-11,
    -0.002_117_147_232_135_5,
    -23.895_741_934_104,
    -5.905_956_432_427e-18,
    -0.000_001_262_180_889_910_1,
    -0.038_946_842_435_739,
    1.125_621_136_045_9e-11,
    -8.231_134_089_799_8,
    1.980_971_280_208_8e-8,
    1.040_696_521_017_4e-19,
    -1.023_474_709_592_9e-13,
    -1.001_817_937_951_1e-9,
    -8.088_290_864_698_5e-11,
    0.106_930_318_794_09,
    -0.336_622_505_741_71,
    8.918_584_535_542_1e-25,
    3.062_931_687_623_2e-13,
    -0.000_004_200_246_769_820_8,
    -5.905_602_968_563_9e-26,
    0.000_003_782_694_761_345_7,
    -1.276_860_893_468_1e-15,
    7.308_761_059_506_1e-29,
    5.541_471_535_077_8e-17,
    -9.436_970_724_121e-7,
];

/// The Region 1 basic equation and its derivatives, up to second order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Gibbs {
    /// `gamma`.
    pub gamma: f64,
    /// `d gamma / d pi`.
    pub gamma_pi: f64,
    /// `d gamma / d tau`.
    pub gamma_tau: f64,
    /// `d^2 gamma / d pi^2`.
    pub gamma_pipi: f64,
    /// `d^2 gamma / d pi d tau`.
    pub gamma_pitau: f64,
    /// `d^2 gamma / d tau^2`.
    pub gamma_tautau: f64,
}

/// The Region 2 ideal part: `gamma0` and its `tau` derivatives. Its `pi` dependence is
/// the single `ln(pi)` term, so its `pi` derivatives are closed forms used inline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ideal {
    /// `gamma0`.
    pub gamma: f64,
    /// `d gamma0 / d tau`.
    pub gamma_tau: f64,
    /// `d^2 gamma0 / d tau^2`.
    pub gamma_tautau: f64,
}

/// The Region 2 residual part and its derivatives, up to second order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Residual {
    /// `gammar`.
    pub gamma: f64,
    /// `d gammar / d pi`.
    pub gamma_pi: f64,
    /// `d gammar / d tau`.
    pub gamma_tau: f64,
    /// `d^2 gammar / d pi^2`.
    pub gamma_pipi: f64,
    /// `d^2 gammar / d pi d tau`.
    pub gamma_pitau: f64,
    /// `d^2 gammar / d tau^2`.
    pub gamma_tautau: f64,
}

/// The Region 1 Gibbs function and its derivatives.
#[must_use]
pub fn gamma1(pi: f64, tau: f64) -> Gibbs {
    let mut g = Gibbs {
        gamma: 0.0,
        gamma_pi: 0.0,
        gamma_tau: 0.0,
        gamma_pipi: 0.0,
        gamma_pitau: 0.0,
        gamma_tautau: 0.0,
    };
    for i in 0..34 {
        let p1 = (7.1 - pi).powi(I1[i]);
        let t1 = (tau - 1.222).powi(J1[i]);
        g.gamma += N1[i] * p1 * t1;
        g.gamma_pi += -N1[i] * I1[i] as f64 * (7.1 - pi).powi(I1[i] - 1) * t1;
        g.gamma_tau += N1[i] * p1 * J1[i] as f64 * (tau - 1.222).powi(J1[i] - 1);
        g.gamma_pipi += N1[i] * I1[i] as f64 * (I1[i] - 1) as f64 * (7.1 - pi).powi(I1[i] - 2) * t1;
        g.gamma_pitau += -N1[i]
            * I1[i] as f64
            * (7.1 - pi).powi(I1[i] - 1)
            * J1[i] as f64
            * (tau - 1.222).powi(J1[i] - 1);
        g.gamma_tautau +=
            N1[i] * p1 * J1[i] as f64 * (J1[i] - 1) as f64 * (tau - 1.222).powi(J1[i] - 2);
    }
    g
}

/// The Region 2 ideal part.
#[must_use]
pub fn gamma0(pi: f64, tau: f64) -> Ideal {
    let mut g = Ideal {
        gamma: pi.ln(),
        gamma_tau: 0.0,
        gamma_tautau: 0.0,
    };
    for i in 0..9 {
        g.gamma += N0[i] * tau.powi(J0[i]);
        g.gamma_tau += N0[i] * J0[i] as f64 * tau.powi(J0[i] - 1);
        g.gamma_tautau += N0[i] * J0[i] as f64 * (J0[i] - 1) as f64 * tau.powi(J0[i] - 2);
    }
    g
}

/// The Region 2 residual part.
#[must_use]
pub fn gammar(pi: f64, tau: f64) -> Residual {
    let mut g = Residual {
        gamma: 0.0,
        gamma_pi: 0.0,
        gamma_tau: 0.0,
        gamma_pipi: 0.0,
        gamma_pitau: 0.0,
        gamma_tautau: 0.0,
    };
    for i in 0..43 {
        g.gamma += NR[i] * pi.powi(IR[i]) * (tau - 0.5).powi(JR[i]);
        g.gamma_pi += NR[i] * IR[i] as f64 * pi.powi(IR[i] - 1) * (tau - 0.5).powi(JR[i]);
        g.gamma_tau += NR[i] * pi.powi(IR[i]) * JR[i] as f64 * (tau - 0.5).powi(JR[i] - 1);
        g.gamma_pipi += NR[i]
            * IR[i] as f64
            * (IR[i] - 1) as f64
            * pi.powi(IR[i] - 2)
            * (tau - 0.5).powi(JR[i]);
        g.gamma_pitau +=
            NR[i] * IR[i] as f64 * pi.powi(IR[i] - 1) * JR[i] as f64 * (tau - 0.5).powi(JR[i] - 1);
        g.gamma_tautau += NR[i]
            * pi.powi(IR[i])
            * JR[i] as f64
            * (JR[i] - 1) as f64
            * (tau - 0.5).powi(JR[i] - 2);
    }
    g
}

/// The saturation pressure in MPa, the Region 4 equation.
#[must_use]
pub fn saturation_pressure(t: f64) -> f64 {
    let theta = t - 0.238_555_575_678_49 / (t - 650.175_348_447_98);
    let a = theta * theta + 1_167.052_145_276_7 * theta - 724_213.167_032_06;
    let b = -17.073_846_940_092 * theta * theta + 12_020.824_702_47 * theta - 3_232_555.032_233_3;
    let c = 14.915_108_613_53 * theta * theta - 4_823.265_736_159_1 * theta + 405_113.405_420_57;
    (2.0 * c / (-b + (b * b - 4.0 * a * c).sqrt())).powi(4)
}

/// The saturation temperature in K, the Region 4 equation.
#[must_use]
pub fn saturation_temperature(p: f64) -> f64 {
    let beta = p.powf(0.25);
    let e = beta * beta - 17.073_846_940_092 * beta + 14.915_108_613_53;
    let f = 1_167.052_145_276_7 * beta * beta + 12_020.824_702_47 * beta - 4_823.265_736_159_1;
    let g = -724_213.167_032_06 * beta * beta - 3_232_555.032_233_3 * beta + 405_113.405_420_57;
    let d = 2.0 * g / (-f - (f * f - 4.0 * e * g).sqrt());
    let s = 650.175_348_447_98 + d;
    (s - (s * s - 4.0 * (-0.238_555_575_678_49 + 650.175_348_447_98 * d)).sqrt()) / 2.0
}

/// The specific thermodynamic properties of one steam state, in IF97's native units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteamProperties {
    /// Compressibility factor, dimensionless.
    pub z: f64,
    /// Specific volume, m³/kg.
    pub v: f64,
    /// Specific internal energy, kJ/kg.
    pub u: f64,
    /// Specific enthalpy, kJ/kg.
    pub h: f64,
    /// Specific entropy, kJ/(kg.K).
    pub s: f64,
    /// Specific isochoric heat capacity, kJ/(kg.K).
    pub cv: f64,
    /// Specific isobaric heat capacity, kJ/(kg.K).
    pub cp: f64,
    /// Specific Gibbs energy, kJ/kg.
    pub g: f64,
    /// Speed of sound, m/s.
    pub w: f64,
}

/// The specific properties at a pressure (MPa) and temperature (K).
#[must_use]
pub fn properties(p_mpa: f64, t: f64) -> SteamProperties {
    let ts = saturation_temperature(p_mpa);
    if t <= ts {
        // Region 1: subcooled liquid.
        let pi = p_mpa / 16.53;
        let tau = 1386.0 / t;
        let g = gamma1(pi, tau);
        let z = pi * g.gamma_pi;
        let h = R * t * tau * g.gamma_tau;
        let s = R * (tau * g.gamma_tau - g.gamma);
        let cp = -R * tau * tau * g.gamma_tautau;
        let cv = cp + R * (g.gamma_pi - tau * g.gamma_pitau).powi(2) / g.gamma_pipi;
        SteamProperties {
            z,
            v: R * t * z / (p_mpa * 1000.0),
            u: h - R * t * z,
            h,
            s,
            cv,
            cp,
            g: R * t * g.gamma,
            w: {
                let num = 1000.0 * R * t * g.gamma_pi * g.gamma_pi;
                let denom = (g.gamma_pi - tau * g.gamma_pitau).powi(2)
                    / (tau * tau * g.gamma_tautau)
                    - g.gamma_pipi;
                (num / denom).sqrt()
            },
        }
    } else {
        // Region 2: superheated vapour.
        let pi = p_mpa;
        let tau = 540.0 / t;
        let id = gamma0(pi, tau);
        let res = gammar(pi, tau);
        let gamma_pi = 1.0 / pi + res.gamma_pi;
        let gamma_pipi = -1.0 / (pi * pi) + res.gamma_pipi;
        let gamma_tau = id.gamma_tau + res.gamma_tau;
        let gamma_tautau = id.gamma_tautau + res.gamma_tautau;
        let gamma = id.gamma + res.gamma;
        let z = pi * gamma_pi;
        let h = R * t * tau * gamma_tau;
        let s = R * (tau * gamma_tau - gamma);
        let cp = -R * tau * tau * gamma_tautau;
        let cv = cp + R * (gamma_pi - tau * res.gamma_pitau).powi(2) / gamma_pipi;
        SteamProperties {
            z,
            v: R * t * z / (p_mpa * 1000.0),
            u: h - R * t * z,
            h,
            s,
            cv,
            cp,
            g: R * t * gamma,
            w: {
                let num = 1000.0
                    * R
                    * t
                    * (1.0 + 2.0 * pi * res.gamma_pi + pi * pi * res.gamma_pi * res.gamma_pi);
                let denom = (1.0 - pi * pi * res.gamma_pipi)
                    + (1.0 + pi * res.gamma_pi - tau * pi * res.gamma_pitau).powi(2)
                        / (tau * tau * gamma_tautau);
                (num / denom).sqrt()
            },
        }
    }
}
