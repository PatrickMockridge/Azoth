//! The projection between a flowsheet and a DEXPI 2.0 Process document.
//!
//! [`export`] writes a flowsheet - and, unless the exchange is topology only, its run - as a
//! Process model. [`import`] reads one back as a scaffold, naming everything it could not
//! recover. Neither direction invents a value: a quantity the run did not compute is omitted,
//! a fluid DEXPI does not carry is left undeclared, and a unit with no reviewed mapping is
//! refused by name.

pub mod export;
pub mod import;

/// The `Azoth/…` property names the extension mode writes.
///
/// **An extension, not part of DEXPI.** A `Data property` is an XSD `nameReference`, not a free
/// string: a bare name may contain no dot, and a dot is legal only after a `prefix/`. The slash
/// form is therefore the schema-valid spelling of a vendor property, and a consumer that does
/// not know these ignores them. They are how azoth's own round trip recovers what the standard
/// model cannot carry - the exact palette unit, the port name, the endpoints, and the recycle.
pub mod extensions {
    /// The palette id of the unit a process step stands for.
    pub const UNIT_OP: &str = "Azoth/UnitOp";
    /// The azoth instance id, where it differs from the step's identifier.
    pub const INSTANCE_ID: &str = "Azoth/InstanceId";
    /// The palette port name a material port stands for.
    pub const PORT: &str = "Azoth/Port";
    /// The index of a `many` outlet's stream.
    pub const INDEX: &str = "Azoth/Index";
    /// The azoth endpoint a stream is produced at.
    pub const FROM: &str = "Azoth/From";
    /// The azoth endpoint a stream is consumed at.
    pub const TO: &str = "Azoth/To";
    /// The tear name of a projected recycle.
    pub const RECYCLE: &str = "Azoth/Recycle";
    /// The flowsheet's own id, which may carry dots the `Model` name may not.
    pub const FLOWSHEET_ID: &str = "Azoth/FlowsheetId";
    /// The flowsheet's display name.
    pub const FLOWSHEET_NAME: &str = "Azoth/FlowsheetName";
}
