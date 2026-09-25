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
  };
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

export interface StreamRecord {
  n: Quantity;
  z: number[];
  P: Quantity;
  T: Quantity;
  h: Quantity;
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
  tears: TearRecord[];
}

/** Everything one call answers with. */
export interface Envelope {
  ok: boolean;
  dirty: boolean;
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
  tools?: Tool[];
}

/** One edit, as the command model reads it. */
export type Command = { command: string } & Record<string, unknown>;
