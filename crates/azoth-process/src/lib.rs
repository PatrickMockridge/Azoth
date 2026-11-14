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
pub mod compressor_chart;
pub mod executor;
pub mod expander_chart;
pub mod flowsheet;
pub mod kernels;
pub mod load;
pub mod middleware;
pub mod model_gen;
pub mod model_inputs_gen;
pub mod model_outputs_gen;
pub mod model_records_gen;
pub mod models;
pub mod order;
pub mod palette_gen;
pub mod pump_chart;
pub mod reactor;
pub mod recycle;
pub mod rotating_design;
pub mod segment;
pub mod separator_design;
pub mod stream;
pub mod surge_curve;
pub mod unit_op;
pub mod unported;
pub mod unported_gen;
pub mod valve;

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
    AbsorptionColumnResult, AntiSurgeRecycleResult, BleveResult, ComponentSplitterResult,
    CompressorResult, CoolerResult, DistillationColumnResult, EjectorResult, EventTreeResult,
    ExpanderResult, FilterResult, FlareResult, GasDispersionResult, GasScrubberResult,
    GaussianPlumeResult, GibbsReactorResult, HeatExchangerResult, HeaterResult, HemReleaseResult,
    JetFireResult, LeakBlowdownResult, ManifoldResult, MapTurboExpanderResult, MixerResult,
    PackedColumnResult, PipeResult, PlugFlowReactorResult, PoolFireResult, ProbitResult,
    PumpResult, RateBasedPackedColumnResult, SeparatorPerformanceResult, SeparatorResult,
    ShortcutDistillationColumnResult, SplitterResult, StirredTankReactorResult,
    StrippingColumnResult, TankResult, ThreePhaseSeparatorResult, ThrottlingValveResult,
    TurboExpanderResult, VceResult, absorption_column, anti_surge_recycle, bleve,
    component_splitter, compressor, cooler, distillation_column, ejector, event_tree, expander,
    filter, flare, gas_dispersion, gas_scrubber, gaussian_plume, gibbs_reactor, heat_exchanger,
    heater, hem_release, jet_fire, leak_blowdown, manifold, map_turbo_expander_compressor, mixer,
    packed_column, pipe, plug_flow_reactor, pool_fire, probit, pump, rate_based_packed_column,
    separator, separator_capacity_assessment, separator_performance, shortcut_distillation_column,
    splitter, stirred_tank_reactor, stripping_column, tank, three_phase_separator,
    throttling_valve, turbo_expander_compressor, vce,
};
pub use stream::Stream;
pub use unit_op::{Param, Source, UnitOpSpec};
