/**
 * What crosses the boundary, as TypeScript.
 *
 * **These are mirrors, and the mirrors are checked rather than trusted.** Every shape below is a
 * Rust struct in `azoth_process::middleware`, and `test/fixtures.test.ts` parses the documents the
 * CLI emitted into these types - so a field renamed in Rust fails a test here instead of
 * producing `undefined` in a widget, which is what the same mistake would otherwise look like.
 *
 * Nothing here re-derives anything the library already decided: the graph is xyflow's own
 * node/edge shape because the projection emits that, and a diagnostic carries a `target` because
 * the checker derives one, so this layer is a reading of the document and not a second opinion
 * about it.
 */

/** Which of the three kinds of node a target or a graph node is. */
export type Role = "instance" | "input" | "product";

export type Severity = "error" | "warning";

/** What to put a mark on, which is not the same question as where in the text. */
export type Target =
  | { kind: "node"; role: Role; id: string }
  | { kind: "handle"; role: Role; node: string; port: string; index: number | null }
  | { kind: "parameter"; node: string; name: string }
  | { kind: "edge"; from: string; to: string }
  | { kind: "endpoint"; endpoint: string }
  | {
      kind: "palette";
      id: string;
      port: string | null;
      parameter: string | null;
      field: string | null;
    }
  | { kind: "document" };

/** Which of the class's two orders the next run takes. */
export type ExecutionOrder = "insertion" | "topological";

/** One diagnostic, as the checker wrote it. */
export interface Diagnostic {
  /** The variant's own name, snake-cased and stable — switch on this, not on `message`. */
  code: string;
  severity: Severity;
  section: string;
  path: string;
  target: Target;
  message: string;
  detail: Record<string, unknown>;
}

export interface Position {
  x: number;
  y: number;
}

/** One port, and the handle ids a stream on it is addressed by. */
export interface Handle {
  name: string;
  multiplicity: "one" | "many";
  handles: string[];
}

export interface Ports {
  inlets: Handle[];
  outlets: Handle[];
}

/** A feed's declared record, in the units the schema declares. */
export interface InputRecord {
  components: string[];
  n: number;
  z: number[];
  P: number;
  T: number;
}

export interface NodeData {
  name: string;
  unit?: string;
  unit_name?: string;
  parameters?: Record<string, unknown>;
  ports?: Ports;
  input?: InputRecord;
}

/** A node, in xyflow's shape — which is the shape the projection emits. */
export interface GraphNode {
  id: string;
  type: "unit_op" | "stream";
  role: Role;
  position: Position;
  data: NodeData;
}

export interface GraphEdge {
  id: string;
  source: string;
  sourceHandle: string;
  target: string;
  targetHandle: string;
  data: {
    kind: "connection" | "recycle";
    from: string;
    to: string;
    /** The session path of the stream this edge carries. */
    path: string;
    /** A tear's seven convergence settings; absent on a connection, which has none. */
    settings?: RecycleSettings;
  };
}

/**
 * A tear's seven settings, as the document states them.
 *
 * **`null` is a silence, not the class's default.** A document that states none of them — the
 * shipped `demo.toml` states none — answers `null` for all seven, and a panel that showed the
 * defaults and wrote them back would turn a silence into a pinned number.
 */
export interface RecycleSettings {
  flow_tolerance: number | null;
  composition_tolerance: number | null;
  temperature_tolerance: number | null;
  pressure_tolerance: number | null;
  max_iterations: number | null;
  minimum_flow: number | null;
  acceleration_method: string | null;
}

export interface Graph {
  id: string;
  name: string;
  nodes: GraphNode[];
  edges: GraphEdge[];
}

/** A scalar with the unit the codec chose for it. */
export interface Quantity {
  magnitude_si: number;
  unit: string;
}

/**
 * A stream, on the palette's five declared fields and the three the run knew.
 *
 * **`mass_flow` and `molar_mass` are `null` where the library cannot weigh the fluid** — a
 * component built from critical constants alone carries no molar mass, and a record that defaulted
 * it to zero would report a mass flow of zero for a real stream. **`vapour_fraction` is `null`
 * where nothing flashed the stream**: it is the one field a record cannot derive, so its absence is
 * a fact about the run and not a zero.
 */
export interface StreamRecord {
  n: Quantity;
  z: number[];
  P: Quantity;
  T: Quantity;
  h: Quantity;
  mass_flow: Quantity | null;
  molar_mass: Quantity | null;
  vapour_fraction: number | null;
}

export interface Residuals {
  flow: number;
  composition: number;
  temperature: number;
  pressure: number;
}

export interface TearRecord {
  stream: string;
  iterations: number;
  solved: boolean;
  active: boolean;
  residuals: Residuals | null;
}

/** The run, as `executor::json` the one codec wrote it. */
export interface SessionReport {
  flowsheet: string;
  converged: boolean;
  iterations: number;
  streams: Record<string, StreamRecord>;
  /**
   * Every unit operation's own answer, by instance id — **the numbers on no outlet stream**.
   *
   * A unit op whose whole answer is its streams has no entry, which is a statement about the
   * arithmetic rather than a gap: a mixer's result *is* its mixed stream. What appears is what a
   * kernel computed *beside* its outlets — a duty, a tray profile, a conversion, a convergence —
   * under the field names that operation's own model declares, so a reader can name a column's
   * `condenser_duty` without a table per unit operation.
   */
  results: Record<string, UnitResult>;
  /** Every declared tear, in declaration order. */
  tears: TearRecord[];
}

/**
 * One unit operation's published result.
 *
 * **The shape is the model's and not this front end's.** The Rust half derives `Serialize` on the
 * struct its case publishes and the Python half publishes the same field names, so an entry here
 * is whatever that operation declares — a scalar, a vector, an array of quantities — and nothing
 * in this file enumerates them. `ResultsSheet` draws what it finds.
 */
export type UnitResult = Record<string, unknown>;

/**
 * One instantiated unit operation, and what backs it.
 *
 * **Every key is always present, `null` rather than missing**, because a reader switching on
 * `verification` should not also have to test for the field's existence — the same rule the
 * diagnostics' own records follow.
 */
export interface UnitProvenance {
  /** The document's own name for this instance, as `Envelope.paths` spells it. */
  instance: string;
  /** The palette entry it was added as, e.g. `unit_ops.pump`. */
  unit: string;
  /** The registered model id, where the palette entry has one. */
  model: string | null;
  /** The provenance line the palette carries, e.g. the NeqSim class. */
  source: string | null;
  /** Whether the executor has a kernel for it, and why not where it does not. */
  runnable: boolean;
  refusal: string | null;
  /**
   * The model's spec and Rust kernel hashes, and how far its answer is checked.
   *
   * **Absent together and absent rather than defaulted**: a palette entry with no model has no
   * calculation to describe, and `null` says so where an empty string would look like a hash
   * nobody could compare.
   */
  spec_sha256: string | null;
  rust_sha256: string | null;
  verification: string | null;
}

/** What a document as a whole rests on: the library, its version, and one entry per instance. */
export interface ProvenanceRecord {
  library: string;
  version: string;
  /**
   * The data authority the run used.
   *
   * The editor's two doors bind the compiled databank and take no overlay, so this is the
   * shipped standard keycard. A reader reproducing a number with its own data needs to know
   * which authority produced this one.
   */
  keycard: string;
  /** One entry per instantiated unit operation, **in the document's own order**. */
  units: UnitProvenance[];
}

/** Everything one call answers with. */
export interface Envelope {
  ok: boolean;
  dirty: boolean;
  /** Which of the class's two orders the next run takes. */
  execution_order: ExecutionOrder;
  flowsheet: {
    id: string;
    name: string;
    document: string;
    graph: Graph;
  };
  diagnostics: Diagnostic[];
  paths: string[];
  session: SessionReport | null;
  run_error: string | null;
  /** What the document rests on, which is a fact about the *run* and not about any stream. */
  provenance: ProvenanceRecord;
}

/** An enum parameter's allowed values live in `values`; this is what chooses a control. */
export type Kind =
  | "quantity"
  | "vector"
  | "boolean"
  | "enum"
  | "string"
  | "components"
  | "matrix"
  | "fitting_list"
  | "record_list"
  | "vector_string"
  | "unknown";

export interface FormField {
  name: string;
  dimension: string;
  shape: "scalar" | "vector";
}

export interface FormPort {
  name: string;
  direction: "in" | "out";
  multiplicity: "one" | "many";
  fields: FormField[];
}

export interface FormRange {
  min: number | null;
  min_inclusive: boolean;
  max: number | null;
  max_inclusive: boolean;
  equals: number | null;
  band: "outside" | "inside";
  severity: Severity;
  code: string;
  rationale: string;
}

/** One parameter, as a field on a form. */
export interface FormParameter {
  name: string;
  kind: Kind;
  required: boolean;
  unit: string | null;
  dimension: string | null;
  values: string[];
  description: string;
  range: FormRange[];
}

/** One unit operation, as a form. */
export interface Form {
  id: string;
  name: string;
  source: string | null;
  model: string | null;
  /**
   * The palette directory the entry was read from, e.g. `two_port`.
   *
   * **The grouping a palette panel draws, and the only place it is stated.** An id is
   * `unit_ops.<leaf>`, so splitting the id takes the leaf for a family and draws one group; and
   * `source` is NeqSim's taxonomy rather than this palette's, where `cooler` is a `two_port` entry
   * inside NeqSim's `heatexchanger/` directory.
   */
  family: string | null;
  runnable: boolean;
  refusal: string | null;
  ports: FormPort[];
  parameters: FormParameter[];
  /**
   * The model's declared outputs, in the spec's order.
   *
   * **The other half of the declaration.** Each carries the unit a bare float is in, the kind,
   * the spec's own sentence and whether it exists only when the design it belongs to was
   * configured. A result sheet explains a value from this; it is not a second schema.
   */
  outputs: CalcOutput[];
  unmodelled_ranges: string[];
}

export interface Tool {
  name: string;
  description: string;
  input_schema: Record<string, unknown>;
}

/** One input of a registered calculation: the same field a unit op's parameter is. */
export interface CalcInput extends FormParameter {
  /** The model's own default, where it declares one; the sheet shows it as the initial value. */
  default: number | boolean | string | null;
}

/** One output, as a result table labels it. */
export interface CalcOutput {
  name: string;
  kind: Kind;
  unit: string | null;
  dimension: string | null;
  values: string[];
  description: string;
  optional: boolean;
}

/** How far the calculation is externally checked, and the hashes that say which code ran. */
export interface CalcProvenance {
  verification: string;
  spec_sha256: string | null;
  rust_sha256: string | null;
  python_sha256: string | null;
  tests_active: number;
  tests_skipped: number;
  validation_cases: number;
}

/** One registered calculation, as the picker and a sheet read it. */
export interface CalcForm {
  id: string;
  name: string;
  /** What the calculation is for, in the spec's own authored words. */
  summary: string;
  namespace: string;
  available: boolean;
  refusal: string | null;
  equation: string | null;
  assumptions: string[];
  references: string[];
  inputs: CalcInput[];
  outputs: CalcOutput[];
  /**
   * The spec's worked example, in the units the spec declares.
   *
   * **Declared units and not SI**, because it is the spec's own text. The sheet converts it with
   * the catalogue's unit factors when a reader asks to load it.
   */
  worked_example: Record<string, number | boolean | string>;
  provenance: CalcProvenance;
}

/** One value a calculation produced: a quantity, or a native enum/flag/string/vector. */
export type SheetValue = number | boolean | string | Quantity | number[];

/** One warning the calculation attached to its answer. */
export interface CalcWarning {
  code: string;
  message: string;
  field: string | null;
}

/** The result of one sheet's last calculation. */
export interface SheetResult {
  outputs: Record<string, SheetValue>;
  warnings: CalcWarning[];
}

/** Why a sheet cannot run, or why it could not: a structural fact, not a kernel refusal. */
export interface SheetDiagnostic {
  code: string;
  message: string;
  field: string | null;
}

/** One sheet in the book. */
export interface Sheet {
  id: string;
  name: string;
  calculation: string;
  inputs: Record<string, number | boolean | string>;
  /** Every required input is set and the calculation is available. */
  ready: boolean;
  missing: string[];
  /** The result is older than the inputs, or there has never been one. */
  stale: boolean;
  result: SheetResult | null;
  error: { message: string; field: string | null } | null;
  diagnostics: SheetDiagnostic[];
}

/** The calculation book, as the middleware owns it. */
export interface SheetBook {
  id: string;
  name: string;
  /** The TOML a save writes; it rides every envelope, so a save is no call. */
  document: string;
  sheets: Sheet[];
}

/** One call's answer from the sheet workspace. */
export interface SheetEnvelope {
  ok: boolean;
  dirty: boolean;
  book: SheetBook;
}

export interface Catalogue {
  unit_ops: Form[];
  /**
   * The standalone calculations, as the sheet editor reads them.
   *
   * Present on the same catalogue as the palette because both are facts about the *library* and a
   * front end has nowhere else to ask for either. A calculation whose inputs or outputs are tables
   * is carried with `available: false` and a refusal rather than omitted.
   */
  calcs?: CalcForm[];
  /**
   * Every canonical unit, keyed by the string a quantity carries: what it measures, and what one
   * of it is worth in SI.
   *
   * **Computed by the library and not tabulated here.** For a scale, `factor` is `si_factor`,
   * which runs the same conversion a calculation runs; for an affine display unit - a degree
   * Celsius, a degree Fahrenheit - `offset` is the unit's own definition, which neither library
   * exposes as data, and `factor` is its scale. `python/tests/test_units_cross_library.py` holds
   * `pint` to each, at one value for a scale and two for a shift. A front end converting with its
   * own table would be a second answer to a question the library has answered, which is the defect
   * the vocabulary's gate exists to refuse.
   *
   * A scale is the same expression with an offset of zero, so a display has one formula:
   * `value / factor - offset`.
   */
  units?: Record<
    string,
    { dimension: string | null; factor: number | null; offset: number | null }
  >;
  /**
   * The named unit sets a reader may switch between, as the vocabulary declares them.
   *
   * One unit per dimension per set, and a dimension no set names is read in the unit the library
   * computed it in - most of the vocabulary has no engineering alternative.
   */
  unit_sets?: { id: string; name: string; units: Record<string, string> }[];
  /**
   * Every substance the compiled databank resolves, with its molar mass where it has one.
   *
   * **A mass column is the databank's to state.** A stream's record carries the *mixture's* molar
   * mass and not each substance's, so a composition table cannot weigh a component without these -
   * and a molecular weight from anywhere else would be a number the run did not use. The mass is
   * kg/mol, the unit the record's own `molar_mass` crosses in, so the two multiply directly.
   */
  components?: { name: string; molar_mass: number | null }[];
  tools?: Tool[];
}

/**
 * One edit, as the command model reads it.
 *
 * Re-exported from `commands.ts` rather than declared twice: that file is where the fifteen are
 * written out, and a second `{ command: string } & Record<string, unknown>` here would let a
 * widget send anything at all and still typecheck.
 */
export type { EditorCommand as Command } from "./commands";

/** What a drawn object would become in azoth. */
export type CoverageStatus =
  | "modelled"
  | "lossy"
  | "ambiguous"
  | "unmapped"
  | "topology"
  | "control"
  | "drawing_only";

/** The exchange has a third severity the checker does not: information, not a loss. */
export type DexpiSeverity = "error" | "warning" | "info";

/** One finding about a DEXPI document or its projection. */
export interface DexpiDiagnostic {
  severity: DexpiSeverity;
  code: string;
  subject: string;
  message: string;
}

export interface DrawingPoint {
  x: number;
  y: number;
}

export interface DrawingColor {
  r: number;
  g: number;
  b: number;
}

export interface DrawingStroke {
  color: DrawingColor;
  width: number;
  dash: string;
}

/** The five primitives a drawing resolves to; `kind` is the wire's tag. */
export type DrawingPrimitive =
  | { kind: "polyline"; id: string | null; points: DrawingPoint[]; stroke: DrawingStroke }
  | {
      kind: "polygon";
      id: string | null;
      points: DrawingPoint[];
      fill: string;
      stroke: DrawingStroke;
    }
  | {
      kind: "ellipse";
      id: string | null;
      center: DrawingPoint;
      radius_x: number;
      radius_y: number;
      rotation: number;
      fill: string;
      stroke: DrawingStroke;
    }
  | {
      kind: "ellipse_arc";
      id: string | null;
      center: DrawingPoint;
      radius_x: number;
      radius_y: number;
      start_angle: number;
      end_angle: number;
      rotation: number;
      stroke: DrawingStroke;
    }
  | {
      kind: "text";
      id: string | null;
      position: DrawingPoint;
      text: string;
      size: number;
      rotation: number;
      alignment: string;
      color: DrawingColor;
    };

export interface DrawingGroup {
  kind: "node" | "edge" | "annotation";
  id: string | null;
  /** The conceptual object the group draws, or an azoth key for a derived schematic. */
  represents: string | null;
  primitives: DrawingPrimitive[];
}

export interface Drawing {
  name: string;
  units: string | null;
  background: DrawingColor;
  bounds: { min_x: number; min_y: number; max_x: number; max_y: number };
  origin: "file" | "derived";
  groups: DrawingGroup[];
}

export interface ObjectCoverage {
  id: string;
  ty: string;
  label: string;
  unit: string | null;
  azoth: string | null;
  status: CoverageStatus;
  explanation: string;
  alternatives: string[];
  imported: boolean;
  exportable: boolean;
  drawn: boolean;
  primitives: number;
}

export interface CoverageCounts {
  semantic: number;
  modelled: number;
  lossy: number;
  ambiguous: number;
  unmapped: number;
  topology: number;
  control: number;
  drawing_only: number;
  drawn: number;
  not_drawn: number;
  groups: number;
  primitives: number;
  layout_positions: number;
}

export interface Scaffold {
  toml: string;
  flowsheet_id: string;
  instances: number;
  connections: number;
  recycles: number;
  instruments: number;
  signals: number;
  layout_positions: number;
  check_errors: number;
  check_warnings: number;
}

/** A drawing, its coverage, and the scaffold the file would become. */
export interface DexpiView {
  format: "dexpi20" | "proteus";
  profile: "process" | "plant" | null;
  model_name: string;
  model_uri: string;
  drawing: Drawing | null;
  objects: ObjectCoverage[];
  counts: CoverageCounts;
  diagnostics: DexpiDiagnostic[];
  scaffold: Scaffold | null;
}

/** The wasm and HTTP doors' answer to `inspect`: a view, or why there is none. */
export interface DexpiViewResponse {
  ok: boolean;
  view?: DexpiView;
  diagnostics?: DexpiDiagnostic[];
}

/** What an export may do, as the Rust `ExchangeOptions` reads it. */
export interface ExchangeOptions {
  allow_lossy: boolean;
  project_recycles: boolean;
  azoth_extensions: boolean;
  topology_only: boolean;
  allow_unconverged: boolean;
  drawing: boolean;
  overrides: Record<string, string>;
}

/** The file, or the diagnostics that refused it. */
export interface DexpiExport {
  ok: boolean;
  xml: string | null;
  diagnostics: DexpiDiagnostic[];
  steps: number;
  streams: number;
  items: number;
  segments: number;
  diagram_groups: number;
  diagram_primitives: number;
}
