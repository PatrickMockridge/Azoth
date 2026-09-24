//! `CatalystBed` - the pellet, the bed and the pressure it costs.
//!
//! A reactor's fixed bed enters the march in three places, and all three are here: the bed's
//! **bulk density** and **specific surface area** convert a rate from a catalyst basis to a
//! volumetric one, its **activity factor** scales the rate, and the **Ergun** equation is the
//! pressure row of the derivative. The Thiele modulus and the effectiveness factor are here
//! too, but they are **off by default** in the reactor
//! (`catalystEffectivenessEnabled = false`), so a caller reaches them only by asking.
//!
//! **The defaults are not neutral.** A default bed is `800 kg/m³` of `3 mm` pellets at a void
//! fraction of `0.40`, which is a real pressure drop; and its activity factor of `1.0` means an
//! unconfigured bed scales a rate by one rather than by zero. Both are the class's own numbers.

/// A fixed bed of catalyst pellets.
///
/// The fields and defaults are `CatalystBed.java:50-71`'s.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CatalystBed {
    /// The pellet diameter the Ergun equation and the Thiele modulus read, m.
    pub particle_diameter: f64,
    /// The bed's void fraction ε.
    pub void_fraction: f64,
    /// The bed's bulk density, kg/m³: what converts a mass-basis rate to a volumetric one.
    pub bulk_density: f64,
    /// The pellet's own density, kg/m³. Read only by the mechanical design.
    pub particle_density: f64,
    /// The pellet's internal porosity εₚ, for the effective diffusivity.
    pub particle_porosity: f64,
    /// The pore tortuosity τ, for the effective diffusivity.
    pub tortuosity: f64,
    /// The specific surface area, m²/kg: what converts an area-basis rate to a volumetric one.
    pub specific_surface_area: f64,
    /// The bed's activity factor, dimensionless in `[0, 1]`.
    pub activity_factor: f64,
}

impl Default for CatalystBed {
    fn default() -> Self {
        Self {
            particle_diameter: 0.003,
            void_fraction: 0.40,
            bulk_density: 800.0,
            particle_density: 1200.0,
            particle_porosity: 0.50,
            tortuosity: 3.0,
            specific_surface_area: 150_000.0,
            activity_factor: 1.0,
        }
    }
}

impl CatalystBed {
    /// The Ergun pressure drop per unit length, Pa/m, positive for a loss.
    ///
    /// `150·μ·(1-ε)²·|u| / (dₚ²·ε³) + 1.75·ρ·(1-ε)·|u|² / (dₚ·ε³)` - the Blake-Kozeny viscous
    /// term and the Burke-Plummer inertial one. **The velocity's sign is discarded**, because
    /// the class writes `Math.abs(superficialVelocity)` and the drop is a loss whichever way
    /// the gas goes.
    #[must_use]
    pub fn pressure_drop(
        &self,
        superficial_velocity: f64,
        gas_density: f64,
        gas_viscosity: f64,
    ) -> f64 {
        let eps = self.void_fraction;
        let dp = self.particle_diameter;
        let u = superficial_velocity.abs();
        let one_minus_eps = 1.0 - eps;
        let eps_cubed = eps * eps * eps;

        let viscous =
            150.0 * gas_viscosity * one_minus_eps * one_minus_eps * u / (dp * dp * eps_cubed);
        let inertial = 1.75 * gas_density * one_minus_eps * u * u / (dp * eps_cubed);
        viscous + inertial
    }

    /// The generalized Thiele modulus for a first-order reaction in a spherical pellet, `φ`.
    ///
    /// `φ = (R/3)·√(|kᵥ| / D_eff)` with `R` the pellet *radius*, and the factor of three the
    /// conversion from a radius-based modulus to the generalized one. **A non-positive
    /// diffusivity gives `1.0e6`** rather than an infinity - the class's own sentinel, which
    /// then lands in the asymptotic effectiveness factor.
    #[must_use]
    pub fn thiele_modulus(&self, volumetric_rate_constant: f64, effective_diffusivity: f64) -> f64 {
        if effective_diffusivity <= 0.0 {
            return 1.0e6;
        }
        let radius = self.particle_diameter / 2.0;
        (radius / 3.0) * (volumetric_rate_constant.abs() / effective_diffusivity).sqrt()
    }

    /// The internal effectiveness factor η for a spherical pellet, `(1/φ)·[1/tanh(3φ) - 1/(3φ)]`.
    ///
    /// Two branches replace the formula at its ends rather than letting it lose precision: `φ`
    /// below `0.01` is `1.0` (no diffusion limitation) and above `500` is `1/(3φ)` (the
    /// asymptotic limit). Between the two it is the class's coth expression.
    #[must_use]
    pub fn effectiveness_factor(&self, thiele_modulus: f64) -> f64 {
        if thiele_modulus < 0.01 {
            return 1.0;
        }
        if thiele_modulus > 500.0 {
            return 1.0 / (3.0 * thiele_modulus);
        }
        let three_phi = 3.0 * thiele_modulus;
        (1.0 / thiele_modulus) * (1.0 / three_phi.tanh() - 1.0 / three_phi)
    }

    /// `D_eff = D·εₚ/τ`, the diffusivity inside the pellet's pores.
    #[must_use]
    pub fn effective_diffusivity(&self, molecular_diffusivity: f64) -> f64 {
        molecular_diffusivity * self.particle_porosity / self.tortuosity
    }
}
