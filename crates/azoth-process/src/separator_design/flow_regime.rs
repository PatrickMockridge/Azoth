//! `MultiphaseFlowRegime`: the inlet pipe's regime map and its droplet distribution.
//!
//! The map is the class's own simplified one: a Mandhane-Gregory-Aziz-style velocity map
//! horizontally, and three velocity criteria vertically (Taylor-bubble, flooding, and the
//! Kutateladze churn-annular transition). The generated distribution is Rosin-Rammler, one
//! correlation per regime, and the entrained-liquid fraction is Oliemans' `tanh` form.
//!
//! The source is `process/equipment/separator/entrainment/MultiphaseFlowRegime.java`.

/// The class's two-phase gas-liquid regimes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowRegime {
    StratifiedSmooth,
    StratifiedWavy,
    Slug,
    Plug,
    Annular,
    AnnularMist,
    DispersedBubble,
    Churn,
    Bubble,
}

impl FlowRegime {
    /// The class's enum spelling.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::StratifiedSmooth => "STRATIFIED_SMOOTH",
            Self::StratifiedWavy => "STRATIFIED_WAVY",
            Self::Slug => "SLUG",
            Self::Plug => "PLUG",
            Self::Annular => "ANNULAR",
            Self::AnnularMist => "ANNULAR_MIST",
            Self::DispersedBubble => "DISPERSED_BUBBLE",
            Self::Churn => "CHURN",
            Self::Bubble => "BUBBLE",
        }
    }
}

/// The inputs the regime map reads.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegimeInput {
    pub gas_density: f64,
    pub liquid_density: f64,
    pub liquid_viscosity: f64,
    pub surface_tension: f64,
    pub pipe_diameter: f64,
    pub gas_superficial_velocity: f64,
    pub liquid_superficial_velocity: f64,
    pub horizontal: bool,
}

/// `predict`: the regime and its Rosin-Rammler `(d63.2, spread)`.
#[must_use]
pub fn predict(input: &RegimeInput) -> (FlowRegime, (f64, f64)) {
    let regime = if input.horizontal {
        horizontal(input)
    } else {
        vertical(input)
    };
    (regime, distribution(regime, input))
}

fn horizontal(input: &RegimeInput) -> FlowRegime {
    let vsg = input.gas_superficial_velocity;
    let vsl = input.liquid_superficial_velocity;

    if vsl < 0.01 {
        if vsg < 3.0 {
            return FlowRegime::StratifiedSmooth;
        } else if vsg < 8.0 {
            return FlowRegime::StratifiedWavy;
        }
        return FlowRegime::Annular;
    }
    if vsl > 3.0 {
        if vsg > 10.0 {
            return FlowRegime::AnnularMist;
        }
        return FlowRegime::DispersedBubble;
    }
    if vsg < 0.3 {
        return FlowRegime::Plug;
    }
    if vsg < 3.0 {
        return FlowRegime::Slug;
    }
    if vsg < 15.0 {
        if vsl < 0.1 {
            return FlowRegime::StratifiedWavy;
        }
        return FlowRegime::Slug;
    }
    if vsl < 0.5 {
        return FlowRegime::Annular;
    }
    FlowRegime::AnnularMist
}

fn vertical(input: &RegimeInput) -> FlowRegime {
    let vsg = input.gas_superficial_velocity;
    let vsl = input.liquid_superficial_velocity;
    let delta = input.liquid_density - input.gas_density;
    let g = 9.81;
    let bubble_slug = 0.25 + 0.35 * (g * input.pipe_diameter * delta / input.liquid_density).sqrt();
    let slug_churn = (g * input.pipe_diameter).sqrt();
    let kutateladze = 3.1
        * (input.surface_tension * g * delta / (input.gas_density * input.gas_density)).powf(0.25);
    let re_mix = input.liquid_density * (vsg + vsl) * input.pipe_diameter / input.liquid_viscosity;

    if vsl > 3.0 && re_mix > 8000.0 {
        return FlowRegime::DispersedBubble;
    }
    if vsg < bubble_slug {
        return FlowRegime::Bubble;
    }
    if vsg < slug_churn {
        return FlowRegime::Slug;
    }
    if vsg < kutateladze {
        return FlowRegime::Churn;
    }
    if vsl < 0.1 {
        return FlowRegime::Annular;
    }
    FlowRegime::AnnularMist
}

fn distribution(regime: FlowRegime, input: &RegimeInput) -> (f64, f64) {
    if input.surface_tension <= 0.0 || input.gas_density <= 0.0 {
        return (100.0e-6, 2.6);
    }
    let d = input.pipe_diameter;
    let vsg = input.gas_superficial_velocity;
    let vsl = input.liquid_superficial_velocity;
    let weber = input.gas_density * vsg * vsg * d / input.surface_tension;
    let weber = weber.max(1.0);
    let re_liquid = input.liquid_density * vsl * d / input.liquid_viscosity;
    match regime {
        FlowRegime::StratifiedSmooth | FlowRegime::StratifiedWavy => {
            let d_max = (weber.powf(-0.6) * d).clamp(50.0e-6, 5.0e-3);
            (d_max / 2.5, 2.0)
        }
        FlowRegime::Annular | FlowRegime::AnnularMist => {
            let re = re_liquid.max(100.0);
            let d32 = (0.069 * weber.powf(-0.6) * re.powf(0.1) * d).clamp(5.0e-6, 1.0e-3);
            (d32 * 1.5, 2.8)
        }
        FlowRegime::Slug => {
            let d_max = (0.725 * weber.powf(-0.6) * d).clamp(30.0e-6, 3.0e-3);
            (d_max / 3.0, 1.8)
        }
        FlowRegime::Churn => {
            let d_max = (0.5 * weber.powf(-0.6) * d).clamp(20.0e-6, 2.0e-3);
            (d_max / 3.0, 2.2)
        }
        FlowRegime::Plug => ((d * 0.05).min(5.0e-3), 2.0),
        FlowRegime::DispersedBubble | FlowRegime::Bubble => {
            let v_mix = vsg + vsl;
            let weber_liquid = input.liquid_density * v_mix * v_mix * d / input.surface_tension;
            let weber_liquid = weber_liquid.max(1.0);
            let d_max = (0.725 * weber_liquid.powf(-0.6) * d).clamp(100.0e-6, 10.0e-3);
            (d_max / 2.0, 2.5)
        }
    }
}

/// `calcEntrainedLiquidFraction`: Oliemans' `tanh(7.25e-7 We^1.25 Re_l^0.25)`.
#[must_use]
pub fn entrained_liquid_fraction(input: &RegimeInput) -> f64 {
    let weber = if input.surface_tension > 0.0 {
        input.gas_density
            * input.gas_superficial_velocity
            * input.gas_superficial_velocity
            * input.pipe_diameter
            / input.surface_tension
    } else {
        0.0
    };
    let re = (input.liquid_density * input.liquid_superficial_velocity * input.pipe_diameter
        / input.liquid_viscosity)
        .max(1.0);
    (7.25e-7 * weber.powf(1.25) * re.powf(0.25)).tanh()
}

#[cfg(test)]
mod tests {
    use super::{FlowRegime, RegimeInput, predict};

    fn base() -> RegimeInput {
        RegimeInput {
            gas_density: 50.0,
            liquid_density: 750.0,
            liquid_viscosity: 2.0e-3,
            surface_tension: 0.025,
            pipe_diameter: 0.2,
            gas_superficial_velocity: 2.0,
            liquid_superficial_velocity: 2.0,
            horizontal: true,
        }
    }

    /// The enhanced capture's two-phase horizontal row predicts `SLUG`, and the row's
    /// generated distribution is the class's slug correlation.
    #[test]
    fn the_captures_horizontal_regime_is_slug() {
        let (regime, (d632, spread)) = predict(&base());
        assert_eq!(regime, FlowRegime::Slug);
        let weber: f64 = 50.0 * 2.0 * 2.0 * 0.2 / 0.025;
        let expected = (0.725 * weber.powf(-0.6) * 0.2).clamp(30.0e-6, 3.0e-3) / 3.0;
        assert!((d632 - expected).abs() / expected < 1.0e-14);
        assert_eq!(spread, 1.8);
    }

    /// A vertical vessel with a low gas rate is the bubble branch.
    #[test]
    fn a_vertical_low_gas_rate_is_bubble() {
        let input = RegimeInput {
            horizontal: false,
            gas_superficial_velocity: 0.01,
            liquid_superficial_velocity: 0.5,
            ..base()
        };
        assert_eq!(predict(&input).0, FlowRegime::Bubble);
    }
}
