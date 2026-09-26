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
  unmodelled_ranges: string[];
}

export interface Tool {
  name: string;
  description: string;
  input_schema: Record<string, unknown>;
}

export interface Catalogue {
  unit_ops: Form[];
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
