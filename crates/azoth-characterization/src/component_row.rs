//! The component row the fluid-rebuild methods read.
//!
//! `PseudoComponentCombiner` lifts a `SystemInterface` into a base-component map and a list of
//! [`PseudoComponentContribution`]s carrying twenty-one numbers each. A fluid is rebuilt from
//! that list, so the row here is the whole of what one component contributes: the names and
//! amounts a blend merges, and the properties a rebuilt pseudo-component averages.
//!
//! **One row type for both halves.** A base component carries its name and its mole amount and
//! leaves every property absent - its constants come from the databank when a mixture is
//! assembled, not from a table a caller retyped. A pseudo-component carries the properties it
//! was characterised with, because no databank has it.
//!
//! Every field is an SI magnitude, which is the shape a `record_list` input crosses as.

/// One component of one fluid.
#[derive(Debug, Clone, PartialEq)]
pub struct FluidComponent {
    /// Which fluid the row belongs to, 0-based. Fluids are ordered by first appearance, and
    /// the row's own order within a fluid is the order the class reads its components in.
    pub fluid: f64,
    /// The component's name. A base component merges by it; a pseudo-component's is the
    /// name it had in its own fluid, and a rebuilt slate renames it.
    pub name: String,
    /// Whether this is a TBP or plus fraction rather than a base component.
    pub pseudo_component: bool,
    /// Its mole amount.
    pub moles: f64,
    /// Its molar mass, in kg/mol. A pseudo-component must state one: it is what the row's
    /// mass is, and what the sorting key falls back to.
    pub molar_mass: Option<f64>,
    /// Its normal liquid density, in kg/m3. A pseudo-component must state one: a group's
    /// volume is the mass over it.
    pub density: Option<f64>,
    /// Its normal boiling point, in K, and its sorting key where it states one.
    pub normal_boiling_point: Option<f64>,
    /// Its critical temperature, in K.
    pub critical_temperature: Option<f64>,
    /// Its critical pressure, in Pa.
    pub critical_pressure: Option<f64>,
    /// Its acentric factor.
    pub acentric_factor: Option<f64>,
    /// Its critical molar volume, in m3/mol.
    pub critical_volume: Option<f64>,
    /// Its Rackett compressibility, and the CPA form beside it. Both drive the Peneloux
    /// volume translation an assembled component is built with.
    pub racket_z: Option<f64>,
    /// See [`Self::racket_z`].
    pub racket_z_cpa: Option<f64>,
    /// Its parachor, in NeqSim's mixed unit.
    pub parachor: Option<f64>,
    /// Its critical viscosity, in Pa*s.
    pub critical_viscosity: Option<f64>,
    /// Its triple-point temperature, in K.
    pub triple_point_temperature: Option<f64>,
    /// Its heat of fusion, in J/mol.
    pub heat_of_fusion: Option<f64>,
    /// Its ideal-gas enthalpy of formation, in J/mol. The blend averages it; no `Component`
    /// carries it, so the assembler drops it.
    pub ideal_gas_enthalpy_of_formation: Option<f64>,
    /// Its ideal-gas heat-capacity coefficients, `cp = a + b T + c T^2 + d T^3`, in
    /// J/(mol*K**n). A pseudo-component must state all four to become a mixture.
    pub cp_a: Option<f64>,
    /// See [`Self::cp_a`].
    pub cp_b: Option<f64>,
    /// See [`Self::cp_a`].
    pub cp_c: Option<f64>,
    /// See [`Self::cp_a`].
    pub cp_d: Option<f64>,
    /// The fitted `m` of the component's attractive term: NeqSim's `addTBPfraction` installs
    /// Soave's `alpha` with this coefficient rather than deriving it from the acentric factor.
    pub attractive_m: Option<f64>,
}

/// One component of a blended or re-cut fluid, as a method returns it.
///
/// The same row less its fluid index: the answer is one fluid, so which of the inputs a
/// component came from is not a property of the answer.
#[derive(Debug, Clone, PartialEq)]
pub struct CharacterizationComponent {
    /// The component's name.
    pub name: String,
    /// Whether this is a pseudo-component.
    pub pseudo_component: bool,
    /// Its mole amount.
    pub moles: f64,
    /// Its molar mass, in kg/mol. Absent for a base component, whose constants are the
    /// databank's rather than the answer's.
    pub molar_mass: Option<f64>,
    /// Its normal liquid density, in kg/m3.
    pub density: Option<f64>,
    /// Its normal boiling point, in K.
    pub normal_boiling_point: Option<f64>,
    /// Its critical temperature, in K.
    pub critical_temperature: Option<f64>,
    /// Its critical pressure, in Pa.
    pub critical_pressure: Option<f64>,
    /// Its acentric factor.
    pub acentric_factor: Option<f64>,
    /// Its critical molar volume, in m3/mol.
    pub critical_volume: Option<f64>,
    /// Its Rackett compressibility.
    pub racket_z: Option<f64>,
    /// Its Rackett compressibility on the CPA basis.
    pub racket_z_cpa: Option<f64>,
    /// Its parachor.
    pub parachor: Option<f64>,
    /// Its critical viscosity, in Pa*s.
    pub critical_viscosity: Option<f64>,
    /// Its triple-point temperature, in K.
    pub triple_point_temperature: Option<f64>,
    /// Its heat of fusion, in J/mol.
    pub heat_of_fusion: Option<f64>,
    /// Its ideal-gas enthalpy of formation, in J/mol.
    pub ideal_gas_enthalpy_of_formation: Option<f64>,
    /// The constant term of its ideal-gas heat capacity, in J/(mol*K).
    pub cp_a: Option<f64>,
    /// The coefficient of `T`, in J/(mol*K**2).
    pub cp_b: Option<f64>,
    /// The coefficient of `T**2`, in J/(mol*K**3).
    pub cp_c: Option<f64>,
    /// The coefficient of `T**3`, in J/(mol*K**4).
    pub cp_d: Option<f64>,
    /// The fitted `m` of its attractive term.
    pub attractive_m: Option<f64>,
}
