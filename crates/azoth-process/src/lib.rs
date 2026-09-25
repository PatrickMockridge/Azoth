//! The process layer: unit operations and flowsheets, on typed, directional
//! channels.
//!
//! A unit operation is a process with inlets and outlets; a channel carries a
//! field record (names → dimensions) with a polarity; conservation is linearity.
//! This crate states the schema — the palette of unit operations and the
//! flowsheet that wires them — and checks that a flowsheet honours the calculus's
//! rules. **The executor that turns a flowsheet into something that runs is
//! [`executor`]**, and a flowsheet is self-contained: its `[[inputs]]` declare the
//! fluid and the state, so it runs with no argument but the document.

pub mod channel;
pub mod check;
pub mod column;
pub mod executor;
pub mod flowsheet;
pub mod kernels;
pub mod load;
pub mod middleware;
pub mod model_gen;
pub mod model_inputs_gen;
pub mod models;
pub mod order;
pub mod palette_gen;
pub mod reactor;
pub mod recycle;
pub mod stream;
pub mod unit_op;

pub use channel::{Direction, FieldType, Multiplicity, Port, Shape};
pub use check::{
    Diagnostic, Location, NodeRole, Severity, Target, split_node_id, validate, validate_palette,
};
pub use flowsheet::{Connection, Flowsheet, Input, Instance, Layout, Recycle};
pub use load::{load_palette, load_palette_text, parse_flowsheet};
pub use order::{ExecutionOrder, execution_order};
pub use recycle::{Acceleration, RecycleSettings, Residuals};
// **The kernels stay under `kernels` and the ids take the flat names**, which is the split
// the Python package makes too. Two call shapes under one name is not a naming problem to
// work around: a kernel takes and returns a `Stream`, which is what a flowsheet's
// connection carries, while a model takes the record field by field so that a case, a
// cross-impl test and a NeqSim capture can address it.
pub use models::{
    AbsorptionColumnResult, ComponentSplitterResult, CompressorResult, CoolerResult,
    DistillationColumnResult, EjectorResult, ExpanderResult, FilterResult, FlareResult,
    GasScrubberResult, GibbsReactorResult, HeatExchangerResult, HeaterResult, ManifoldResult,
    MixerResult, PackedColumnResult, PipeResult, PlugFlowReactorResult, PumpResult,
    SeparatorResult, ShortcutDistillationColumnResult, SplitterResult, StirredTankReactorResult,
    StrippingColumnResult, TankResult, ThreePhaseSeparatorResult, ThrottlingValveResult,
    absorption_column, component_splitter, compressor, cooler, distillation_column, ejector,
    expander, filter, flare, gas_scrubber, gibbs_reactor, heat_exchanger, heater, manifold, mixer,
    packed_column, pipe, plug_flow_reactor, pump, separator, shortcut_distillation_column,
    splitter, stirred_tank_reactor, stripping_column, tank, three_phase_separator,
    throttling_valve,
};
pub use stream::Stream;
pub use unit_op::{Param, Source, UnitOpSpec};
