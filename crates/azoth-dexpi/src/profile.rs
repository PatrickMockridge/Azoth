//! The DEXPI information-model profile this crate exchanges.
//!
//! Tranche one is the **Process** profile - a PFD/BFD: process steps, material ports, streams
//! and physical state quantities. The Plant profile (P&ID) and the Proteus serialisation are
//! later tranches, and no enum variant is declared for them until they land.

/// The Core information model, at the release this crate pins.
pub const CORE_MODEL: &str = "https://data.dexpi.org/models/2.0.0/Core.xml";

/// The Process information model, at the release this crate pins.
pub const PROCESS_MODEL: &str = "https://data.dexpi.org/models/2.0.0/Process.xml";

/// The profile an exchange is written to or read as.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Profile {
    /// Core + Process: PFD/BFD content.
    Process,
}

impl Profile {
    /// The information models this profile imports, in the order its imports are written.
    #[must_use]
    pub fn imports(self) -> [(&'static str, &'static str); 2] {
        match self {
            Self::Process => [("Core", CORE_MODEL), ("Process", PROCESS_MODEL)],
        }
    }
}

/// The object types the Process profile supports.
///
/// **The reviewed set, not the DEXPI model's whole class hierarchy.** This mirrors NeqSim's
/// `Dexpi20SemanticValidator.SUPPORTED_TYPES`, restricted to the Core and Process entries: a
/// type outside it has not been reviewed for this port, and the validator refuses it rather
/// than guessing it will be understood.
pub const SUPPORTED_TYPES: &[&str] = &[
    "Core/EngineeringModel",
    "Core/QualifiedValue",
    "Core/PhysicalQuantities.PhysicalQuantity",
    "Process/ProcessModel",
    "Process/Process.MaterialPort",
    "Process/Process.Stream",
    "Process/Process.Source",
    "Process/Process.Sink",
    "Process/Process.TransportingFluids",
    "Process/Process.Compressing",
    "Process/Process.Pumping",
    "Process/Process.Distilling",
    "Process/Process.SeparatingByGravity",
    "Process/Process.Cooling",
    "Process/Process.ExchangingThermalEnergy",
    "Process/Process.HeatingInFurnace",
    "Process/Process.StoringFluids",
    "Process/Process.RegulatingFlow",
    "Process/Process.MixingSimple",
    "Process/Process.SplittingMaterial",
];

/// The port directions a `MaterialPort` may carry.
pub const DIRECTION_INLET: &str = "Process/Enumerations.PortDirection.Inlet";
/// The outlet direction.
pub const DIRECTION_OUTLET: &str = "Process/Enumerations.PortDirection.Outlet";

/// Whether a type is in the reviewed Process profile.
#[must_use]
pub fn supports(ty: &str) -> bool {
    SUPPORTED_TYPES.contains(&ty)
}
