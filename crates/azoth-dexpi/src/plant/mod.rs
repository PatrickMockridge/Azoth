//! The projection between a flowsheet and a DEXPI 2.0 **Plant** exchange: the P&ID.
//!
//! [`export`] writes the tagged items and the piping between them; [`import`] reads a file back
//! as a flowsheet scaffold. The Plant model has no process state - no mass flow, pressure or
//! temperature - so neither direction carries quantities. What it carries is the part a PFD
//! loses: which line a stream belongs to, which nozzle it leaves by, and the tag of the item it
//! leaves.
//!
//! **A Plant exchange is not a less-detailed Process one.** They are two official information
//! models, and this module is the second projection beside [`crate::process`].

pub mod export;
pub mod import;

use crate::mapping::Disposition;
use crate::plant_mapping_gen::{PLANT_CLASS_CANDIDATES, PLANT_UNITS};

/// One palette entry's Plant mapping - the write side and the read side, which are not the
/// same set.
#[derive(Debug, Clone, Copy)]
pub struct PlantMapping {
    pub azoth_id: &'static str,
    /// What writing the unit as a tagged item costs.
    pub write: Disposition,
    /// The native `Plant/ProcessEquipment.*` class the writer emits; empty for a refusal.
    pub plant_types: &'static [&'static str],
    /// The bare class names a reader accepts, native or Proteus.
    pub reader_types: &'static [&'static str],
    /// What a lossy write loses, or why a write is refused.
    pub note: &'static str,
}

/// Every Plant mapping row, in palette-id order.
#[must_use]
pub fn units() -> &'static [PlantMapping] {
    PLANT_UNITS
}

/// The row for one azoth palette id.
#[must_use]
pub fn for_azoth(id: &str) -> Option<&'static PlantMapping> {
    PLANT_UNITS.iter().find(|row| row.azoth_id == id)
}

/// The azoth unit a class name stands for. A native type is spelled
/// `Plant/ProcessEquipment.<name>`, a Proteus one `<name>`; the class is the last segment.
#[must_use]
pub fn for_class(class: &str) -> Option<&'static str> {
    PLANT_CLASS_CANDIDATES
        .iter()
        .find_map(|(name, id)| (*name == class).then_some(*id))
}

/// The class name of a DEXPI type: the last dotted segment. `Plant/Piping.GlobeValve` is
/// `GlobeValve`; a bare Proteus `ComponentClass` is itself.
#[must_use]
pub fn class_of(ty: &str) -> &str {
    ty.rsplit('.').next().unwrap_or(ty)
}

/// The Plant mapping's counts: faithful, lossy and refused writes, and readable rows.
#[must_use]
pub fn counts() -> (usize, usize, usize, usize) {
    (
        crate::plant_mapping_gen::PLANT_FAITHFUL_COUNT,
        crate::plant_mapping_gen::PLANT_LOSSY_COUNT,
        crate::plant_mapping_gen::PLANT_REFUSED_COUNT,
        crate::plant_mapping_gen::PLANT_READABLE_COUNT,
    )
}

/// The line label a flowsheet with no `line` fields is written under.
///
/// **It is azoth's own placeholder, and import recognises it.** DEXPI wants a
/// `PipingNetworkSystem` per line; a flowsheet that states no lines still has material edges,
/// so they go under this one label, and a reader turns it back into "no line" rather than
/// inventing a line the user never wrote.
pub const DEFAULT_LINE: &str = "AzothProcessConnections";

/// The classes a reader collapses rather than turning into a unit.
///
/// **A pass-through is not a unit.** A reducer, a blind flange, a strainer or an orifice plate
/// is a fitting in a line; the stream crosses it unchanged as far as this model is concerned.
/// The list is reviewed against NeqSim's piping-component table, where these names exist and
/// map to its `Stream` - the same statement this list makes.
pub const PASS_THROUGH_CLASSES: &[&str] = &[
    "Pipe",
    "PipeReducer",
    "BlindFlange",
    "Strainer",
    "FlameArrestor",
    "OrificePlate",
    "FlowMeter",
    "VenturiMeter",
];
