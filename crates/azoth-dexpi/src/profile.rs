//! The DEXPI information-model profiles this crate exchanges: the P&ID and the PFD/BFD.
//!
//! A **Process** exchange carries process steps, material ports, streams and physical state
//! quantities; a **Plant** exchange carries tagged plant items, nozzles, piping network
//! systems and segments. They are two official information models of one release, not a
//! detailed and a simple version of one diagram.

use crate::model::Document;

/// The Core information model, at the release this crate pins.
pub const CORE_MODEL: &str = "https://data.dexpi.org/models/2.0.0/Core.xml";

/// The Process information model, at the release this crate pins.
pub const PROCESS_MODEL: &str = "https://data.dexpi.org/models/2.0.0/Process.xml";

/// The Plant information model, at the release this crate pins.
pub const PLANT_MODEL: &str = "https://data.dexpi.org/models/2.0.0/Plant.xml";

/// The profile an exchange is written to or read as.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Profile {
    /// Core + Process: PFD/BFD content.
    Process,
    /// Core + Plant: P&ID content.
    Plant,
}

impl Profile {
    /// The information models this profile imports, in the order its imports are written.
    #[must_use]
    pub fn imports(self) -> [(&'static str, &'static str); 2] {
        match self {
            Self::Process => [("Core", CORE_MODEL), ("Process", PROCESS_MODEL)],
            Self::Plant => [("Core", CORE_MODEL), ("Plant", PLANT_MODEL)],
        }
    }

    /// The profile a document declares through its imports.
    ///
    /// A document importing both is refused: the two models are alternatives, and one file
    /// claiming both has not said which exchange it is.
    #[must_use]
    pub fn of(document: &Document) -> Option<Self> {
        let process = document
            .model
            .imports
            .iter()
            .any(|import| import.source == PROCESS_MODEL);
        let plant = document
            .model
            .imports
            .iter()
            .any(|import| import.source == PLANT_MODEL);
        match (plant, process) {
            (true, false) => Some(Self::Plant),
            (false, true) => Some(Self::Process),
            _ => None,
        }
    }

    /// The object types this profile supports.
    #[must_use]
    pub fn supported_types(self) -> &'static [&'static str] {
        match self {
            Self::Process => SUPPORTED_TYPES,
            Self::Plant => PLANT_SUPPORTED_TYPES,
        }
    }

    /// Whether a type is in this profile's reviewed set.
    ///
    /// **A Plant projection reads more classes than the native writer reviewed.** The writer's
    /// set is narrow because it was reviewed for writing; the reader accepts every class in
    /// `specs/dexpi/plant_units.toml`, and a projection of those classes spells them as native
    /// types. Both lists are this port's reviewed sets, and the validator holds a document to
    /// the union.
    #[must_use]
    pub fn supports(self, ty: &str) -> bool {
        match self {
            Self::Process => SUPPORTED_TYPES.contains(&ty),
            Self::Plant => {
                PLANT_SUPPORTED_TYPES.contains(&ty)
                    || crate::plant_mapping_gen::PLANT_CLASS_TYPES.contains(&ty)
            }
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

/// The object types the Plant profile supports.
///
/// The Plant half of NeqSim's reviewed set. Instrumentation types are readable because a real
/// P&ID carries them - this port does not yet model instruments, and the importer reports an
/// instrument as unplaced rather than the validator refusing the file.
pub const PLANT_SUPPORTED_TYPES: &[&str] = &[
    "Core/EngineeringModel",
    "Core/QualifiedValue",
    "Core/PhysicalQuantities.PhysicalQuantity",
    "Plant/PlantModel",
    "Plant/ProcessEquipment.ProcessEquipment",
    "Plant/ProcessEquipment.CentrifugalCompressor",
    "Plant/ProcessEquipment.CentrifugalPump",
    "Plant/ProcessEquipment.Separator",
    "Plant/ProcessEquipment.AirCoolingSystem",
    "Plant/ProcessEquipment.TubularHeatExchanger",
    "Plant/ProcessEquipment.FiredHeater",
    "Plant/ProcessEquipment.Tank",
    "Plant/ProcessEquipment.Nozzle",
    "Plant/Piping.PipingNetworkSystem",
    "Plant/Piping.PipingNetworkSegment",
    "Plant/Piping.PipingNode",
    "Plant/Piping.Pipe",
    "Plant/Piping.GlobeValve",
    "Plant/Piping.SpringLoadedGlobeSafetyValve",
    "Plant/Piping.SwingCheckValve",
    "Plant/Piping.FlowInPipeOffPageConnector",
    "Plant/Piping.FlowOutPipeOffPageConnector",
    "Plant/Instrumentation.ProcessInstrumentationFunction",
    "Plant/Instrumentation.ProcessSignalGeneratingFunction",
    "Plant/Instrumentation.MeasuringLineFunction",
    "Plant/Instrumentation.SignalConveyingFunction",
    "Plant/Instrumentation.InstrumentationLoopFunction",
    "Plant/Diagram.PlantMetaData",
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
