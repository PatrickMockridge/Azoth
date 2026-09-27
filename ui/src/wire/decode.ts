/**
 * The boundary as a parse.
 *
 * **Every document that crosses the wire is decoded, and none is cast.** `client.ts` and `http.ts`
 * hand this module what the library wrote — a `JSON.parse` result, or the wasm module's own
 * string — and take back a value of the type the rest of the editor dereferences. A cast asserts
 * the shape; a parse *finds out*, and when it is wrong it says which path was wrong and what it
 * wanted there. The difference is what a field renamed in Rust looks like: a cast produces
 * `undefined` in a widget three files away, and a parse names the field at the boundary it crossed.
 *
 * **This module is the other half of the mirror.** `types.ts` is hand-written from the Rust structs
 * and this file is what `test/fixtures.test.ts` checks it against: the CLI's own documents are
 * decoded by the same code a browser runs, so a shape the declaration stopped matching fails the
 * test rather than the reader.
 *
 * **Where it stops, and why.** It decodes every shape a panel *addresses* — the envelope's top
 * level, the graph and its handles, a diagnostic and the target it marks, a stream record and its
 * quantities, a form with its ports and parameters. It does not descend into the four things the
 * library declares open: a diagnostic's `detail`, a node's `parameters`, a tool's `input_schema`
 * and a unit op's published result are `Record<string, unknown>` by declaration, and a shape fixed
 * here would be a second opinion about them rather than a reading of one.
 *
 * ui-ok-file: this module is the decoder the cast rule exempts — a cast inside a checked function is
 * the whole technique of parsing an untyped document, and the exemption is this file and no other.
 */

import type {
  Catalogue,
  Diagnostic,
  Envelope,
  ExecutionOrder,
  Form,
  FormField,
  FormParameter,
  FormPort,
  FormRange,
  Graph,
  GraphEdge,
  GraphNode,
  Handle,
  InputRecord,
  Kind,
  NodeData,
  Ports,
  Position,
  Quantity,
  RecycleSettings,
  Residuals,
  Role,
  SessionReport,
  Severity,
  StreamRecord,
  Target,
  TearRecord,
} from "./types";

/** What crossed was not the shape the declaration says, and this says which path was wrong. */
export class WireError extends Error {
  constructor(
    readonly path: string,
    readonly wanted: string,
    readonly got: unknown,
  ) {
    super(`${path}: wanted ${wanted}, and the library wrote ${describe(got)}`);
    this.name = "WireError";
  }
}

function describe(value: unknown): string {
  if (value === undefined) {
    // The shape a *renamed* field takes, and the one worth saying plainly: the field is not there.
    return "nothing at all";
  }
  if (value === null) {
    return "null";
  }
  if (Array.isArray(value)) {
    return "an array";
  }
  if (typeof value === "object") {
    return "an object";
  }
  return `${typeof value} ${JSON.stringify(value)}`;
}

function fail(path: string, wanted: string, got: unknown): never {
  throw new WireError(path, wanted, got);
}

function obj(value: unknown, path: string): Record<string, unknown> {
  if (typeof value !== "object" || value === null || Array.isArray(value)) {
    fail(path, "an object", value);
  }
  return value as Record<string, unknown>;
}

function str(value: unknown, path: string): string {
  return typeof value === "string" ? value : fail(path, "a string", value);
}

function num(value: unknown, path: string): number {
  return typeof value === "number" ? value : fail(path, "a number", value);
}

function bool(value: unknown, path: string): boolean {
  return typeof value === "boolean" ? value : fail(path, "a boolean", value);
}

function list<T>(value: unknown, path: string, item: (v: unknown, p: string) => T): T[] {
  if (!Array.isArray(value)) {
    fail(path, "an array", value);
  }
  return value.map((entry, index) => item(entry, `${path}[${index}]`));
}

function nullable<T>(
  value: unknown,
  path: string,
  some: (v: unknown, p: string) => T,
): T | null {
  return value === null || value === undefined ? null : some(value, path);
}

/** One of a closed set of strings, which is how an enum crosses. */
function oneOf<T extends string>(value: unknown, path: string, allowed: readonly T[]): T {
  const found = str(value, path);
  const match = allowed.find((entry) => entry === found);
  return match ?? fail(path, `one of ${allowed.join(", ")}`, value);
}

/** A map whose keys are the document's, and whose values are a shape. */
function mapOf<T>(
  value: unknown,
  path: string,
  each: (v: unknown, p: string) => T,
): Record<string, T> {
  const entries = Object.entries(obj(value, path));
  const built: Record<string, T> = {};
  for (const [key, entry] of entries) {
    built[key] = each(entry, `${path}.${key}`);
  }
  return built;
}

const ROLES = ["instance", "input", "product"] as const;
const SEVERITIES = ["error", "warning"] as const;
const ORDERS = ["insertion", "topological"] as const;
const MULTIPLICITIES = ["one", "many"] as const;
const KINDS = [
  "quantity",
  "vector",
  "boolean",
  "enum",
  "string",
  "components",
  "matrix",
  "unknown",
] as const;

const role = (value: unknown, path: string): Role => oneOf(value, path, ROLES);
const severity = (value: unknown, path: string): Severity => oneOf(value, path, SEVERITIES);

/**
 * The role a node id names, or `null` when it names none of the three.
 *
 * **The parse `nodes.ts` used to cast.** `{role}:{name}` is a projection's convention and not a
 * guarantee about the string, so a reader that asserted the role would invent a command for a node
 * the projection did not make. Returning `null` is the same answer the switch's `default` gives,
 * arrived at by reading rather than by asserting.
 */
export function asRole(value: string): Role | null {
  return ROLES.find((entry) => entry === value) ?? null;
}

function position(value: unknown, path: string): Position {
  const raw = obj(value, path);
  return { x: num(raw.x, `${path}.x`), y: num(raw.y, `${path}.y`) };
}

function quantity(value: unknown, path: string): Quantity {
  const raw = obj(value, path);
  return {
    magnitude_si: num(raw.magnitude_si, `${path}.magnitude_si`),
    unit: str(raw.unit, `${path}.unit`),
  };
}

function handle(value: unknown, path: string): Handle {
  const raw = obj(value, path);
  return {
    name: str(raw.name, `${path}.name`),
    multiplicity: oneOf(raw.multiplicity, `${path}.multiplicity`, MULTIPLICITIES),
    handles: list(raw.handles, `${path}.handles`, str),
  };
}

function ports(value: unknown, path: string): Ports {
  const raw = obj(value, path);
  return {
    inlets: list(raw.inlets, `${path}.inlets`, handle),
    outlets: list(raw.outlets, `${path}.outlets`, handle),
  };
}

function inputRecord(value: unknown, path: string): InputRecord {
  const raw = obj(value, path);
  return {
    components: list(raw.components, `${path}.components`, str),
    n: num(raw.n, `${path}.n`),
    z: list(raw.z, `${path}.z`, num),
    P: num(raw.P, `${path}.P`),
    T: num(raw.T, `${path}.T`),
  };
}

function nodeData(value: unknown, path: string): NodeData {
  const raw = obj(value, path);
  const data: NodeData = { name: str(raw.name, `${path}.name`) };
  if (raw.unit !== undefined) {
    data.unit = str(raw.unit, `${path}.unit`);
  }
  if (raw.unit_name !== undefined) {
    data.unit_name = str(raw.unit_name, `${path}.unit_name`);
  }
  if (raw.ports !== undefined) {
    data.ports = ports(raw.ports, `${path}.ports`);
  }
  if (raw.input !== undefined) {
    data.input = inputRecord(raw.input, `${path}.input`);
  }
  // `parameters` is the library's to shape and this layer's to carry: a document states whatever
  // its own spec declares, so the values are objects by check and opaque by declaration.
  if (raw.parameters !== undefined) {
    data.parameters = obj(raw.parameters, `${path}.parameters`);
  }
  return data;
}

function graphNode(value: unknown, path: string): GraphNode {
  const raw = obj(value, path);
  return {
    id: str(raw.id, `${path}.id`),
    type: oneOf(raw.type, `${path}.type`, ["unit_op", "stream"] as const),
    role: role(raw.role, `${path}.role`),
    position: position(raw.position, `${path}.position`),
    data: nodeData(raw.data, `${path}.data`),
  };
}

function recycleSettings(value: unknown, path: string): RecycleSettings {
  const raw = obj(value, path);
  const measured = [
    "flow_tolerance",
    "composition_tolerance",
    "temperature_tolerance",
    "pressure_tolerance",
    "max_iterations",
    "minimum_flow",
  ] as const satisfies readonly (keyof RecycleSettings)[];
  const settings = {
    acceleration_method: nullable(raw.acceleration_method, `${path}.acceleration_method`, str),
  } as RecycleSettings;
  for (const field of measured) {
    settings[field] = nullable(raw[field], `${path}.${field}`, num);
  }
  return settings;
}

function graphEdge(value: unknown, path: string): GraphEdge {
  const raw = obj(value, path);
  const data = obj(raw.data, `${path}.data`);
  const edge: GraphEdge = {
    id: str(raw.id, `${path}.id`),
    source: str(raw.source, `${path}.source`),
    sourceHandle: str(raw.sourceHandle, `${path}.sourceHandle`),
    target: str(raw.target, `${path}.target`),
    targetHandle: str(raw.targetHandle, `${path}.targetHandle`),
    data: {
      kind: oneOf(data.kind, `${path}.data.kind`, ["connection", "recycle"] as const),
      from: str(data.from, `${path}.data.from`),
      to: str(data.to, `${path}.data.to`),
      path: str(data.path, `${path}.data.path`),
    },
  };
  if (data.settings !== undefined) {
    edge.data.settings = recycleSettings(data.settings, `${path}.data.settings`);
  }
  return edge;
}

function graph(value: unknown, path: string): Graph {
  const raw = obj(value, path);
  return {
    id: str(raw.id, `${path}.id`),
    name: str(raw.name, `${path}.name`),
    nodes: list(raw.nodes, `${path}.nodes`, graphNode),
    edges: list(raw.edges, `${path}.edges`, graphEdge),
  };
}

function target(value: unknown, path: string): Target {
  const raw = obj(value, path);
  const kind = str(raw.kind, `${path}.kind`);
  switch (kind) {
    case "node":
      return { kind, role: role(raw.role, `${path}.role`), id: str(raw.id, `${path}.id`) };
    case "handle":
      return {
        kind,
        role: role(raw.role, `${path}.role`),
        node: str(raw.node, `${path}.node`),
        port: str(raw.port, `${path}.port`),
        index: nullable(raw.index, `${path}.index`, num),
      };
    case "parameter":
      return {
        kind,
        node: str(raw.node, `${path}.node`),
        name: str(raw.name, `${path}.name`),
      };
    case "edge":
      return { kind, from: str(raw.from, `${path}.from`), to: str(raw.to, `${path}.to`) };
    case "endpoint":
      return { kind, endpoint: str(raw.endpoint, `${path}.endpoint`) };
    case "palette":
      return {
        kind,
        id: str(raw.id, `${path}.id`),
        port: nullable(raw.port, `${path}.port`, str),
        parameter: nullable(raw.parameter, `${path}.parameter`, str),
        field: nullable(raw.field, `${path}.field`, str),
      };
    case "document":
      return { kind };
    default:
      return fail(path, "a target kind a canvas can mark", value);
  }
}

function diagnostic(value: unknown, path: string): Diagnostic {
  const raw = obj(value, path);
  return {
    code: str(raw.code, `${path}.code`),
    severity: severity(raw.severity, `${path}.severity`),
    section: str(raw.section, `${path}.section`),
    path: str(raw.path, `${path}.path`),
    target: target(raw.target, `${path}.target`),
    message: str(raw.message, `${path}.message`),
    detail: obj(raw.detail, `${path}.detail`),
  };
}

function streamRecord(value: unknown, path: string): StreamRecord {
  const raw = obj(value, path);
  return {
    n: quantity(raw.n, `${path}.n`),
    z: list(raw.z, `${path}.z`, num),
    P: quantity(raw.P, `${path}.P`),
    T: quantity(raw.T, `${path}.T`),
    h: quantity(raw.h, `${path}.h`),
    mass_flow: nullable(raw.mass_flow, `${path}.mass_flow`, quantity),
    molar_mass: nullable(raw.molar_mass, `${path}.molar_mass`, quantity),
    vapour_fraction: nullable(raw.vapour_fraction, `${path}.vapour_fraction`, num),
  };
}

function residuals(value: unknown, path: string): Residuals {
  const raw = obj(value, path);
  return {
    flow: num(raw.flow, `${path}.flow`),
    composition: num(raw.composition, `${path}.composition`),
    temperature: num(raw.temperature, `${path}.temperature`),
    pressure: num(raw.pressure, `${path}.pressure`),
  };
}

function tearRecord(value: unknown, path: string): TearRecord {
  const raw = obj(value, path);
  return {
    stream: str(raw.stream, `${path}.stream`),
    iterations: num(raw.iterations, `${path}.iterations`),
    solved: bool(raw.solved, `${path}.solved`),
    active: bool(raw.active, `${path}.active`),
    residuals: nullable(raw.residuals, `${path}.residuals`, residuals),
  };
}

function sessionReport(value: unknown, path: string): SessionReport {
  const raw = obj(value, path);
  return {
    flowsheet: str(raw.flowsheet, `${path}.flowsheet`),
    converged: bool(raw.converged, `${path}.converged`),
    iterations: num(raw.iterations, `${path}.iterations`),
    streams: mapOf(raw.streams, `${path}.streams`, streamRecord),
    // **The one shape this file refuses to fix**, because the library does not: a published result
    // is whatever the operation's own model declares, and enumerating it here would be a table per
    // unit operation on the wrong side of the boundary.
    results: mapOf(raw.results, `${path}.results`, (entry, at) => obj(entry, at)),
    tears: list(raw.tears, `${path}.tears`, tearRecord),
  };
}

/** Everything one call answers with. The parse `client.ts` and `http.ts` both end at. */
export function decodeEnvelope(value: unknown): Envelope {
  const path = "the envelope";
  const raw = obj(value, path);
  const flowsheet = obj(raw.flowsheet, `${path}.flowsheet`);
  return {
    ok: bool(raw.ok, `${path}.ok`),
    dirty: bool(raw.dirty, `${path}.dirty`),
    execution_order: oneOf<ExecutionOrder>(
      raw.execution_order,
      `${path}.execution_order`,
      ORDERS,
    ),
    flowsheet: {
      id: str(flowsheet.id, `${path}.flowsheet.id`),
      name: str(flowsheet.name, `${path}.flowsheet.name`),
      document: str(flowsheet.document, `${path}.flowsheet.document`),
      graph: graph(flowsheet.graph, `${path}.flowsheet.graph`),
    },
    diagnostics: list(raw.diagnostics, `${path}.diagnostics`, diagnostic),
    paths: list(raw.paths, `${path}.paths`, str),
    session: nullable(raw.session, `${path}.session`, sessionReport),
    run_error: nullable(raw.run_error, `${path}.run_error`, str),
  };
}

function formField(value: unknown, path: string): FormField {
  const raw = obj(value, path);
  return {
    name: str(raw.name, `${path}.name`),
    dimension: str(raw.dimension, `${path}.dimension`),
    shape: oneOf(raw.shape, `${path}.shape`, ["scalar", "vector"] as const),
  };
}

function formPort(value: unknown, path: string): FormPort {
  const raw = obj(value, path);
  return {
    name: str(raw.name, `${path}.name`),
    direction: oneOf(raw.direction, `${path}.direction`, ["in", "out"] as const),
    multiplicity: oneOf(raw.multiplicity, `${path}.multiplicity`, MULTIPLICITIES),
    fields: list(raw.fields, `${path}.fields`, formField),
  };
}

function formRange(value: unknown, path: string): FormRange {
  const raw = obj(value, path);
  return {
    min: nullable(raw.min, `${path}.min`, num),
    min_inclusive: bool(raw.min_inclusive, `${path}.min_inclusive`),
    max: nullable(raw.max, `${path}.max`, num),
    max_inclusive: bool(raw.max_inclusive, `${path}.max_inclusive`),
    equals: nullable(raw.equals, `${path}.equals`, num),
    band: oneOf(raw.band, `${path}.band`, ["outside", "inside"] as const),
    severity: severity(raw.severity, `${path}.severity`),
    code: str(raw.code, `${path}.code`),
    rationale: str(raw.rationale, `${path}.rationale`),
  };
}

function formParameter(value: unknown, path: string): FormParameter {
  const raw = obj(value, path);
  return {
    name: str(raw.name, `${path}.name`),
    kind: oneOf<Kind>(raw.kind, `${path}.kind`, KINDS),
    required: bool(raw.required, `${path}.required`),
    unit: nullable(raw.unit, `${path}.unit`, str),
    dimension: nullable(raw.dimension, `${path}.dimension`, str),
    values: list(raw.values, `${path}.values`, str),
    description: str(raw.description, `${path}.description`),
    range: list(raw.range, `${path}.range`, formRange),
  };
}

function form(value: unknown, path: string): Form {
  const raw = obj(value, path);
  return {
    id: str(raw.id, `${path}.id`),
    name: str(raw.name, `${path}.name`),
    source: nullable(raw.source, `${path}.source`, str),
    model: nullable(raw.model, `${path}.model`, str),
    family: nullable(raw.family, `${path}.family`, str),
    runnable: bool(raw.runnable, `${path}.runnable`),
    refusal: nullable(raw.refusal, `${path}.refusal`, str),
    ports: list(raw.ports, `${path}.ports`, formPort),
    parameters: list(raw.parameters, `${path}.parameters`, formParameter),
    unmodelled_ranges: list(raw.unmodelled_ranges, `${path}.unmodelled_ranges`, str),
  };
}

/** The palette, the units and the components: what a front end asks the *library* for. */
export function decodeCatalogue(value: unknown): Catalogue {
  const path = "the catalogue";
  const raw = obj(value, path);
  const units = nullable(raw.units, `${path}.units`, (entry, at) =>
    mapOf(entry, at, (unit, unitPath) => {
      const declaration = obj(unit, unitPath);
      return {
        dimension: nullable(declaration.dimension, `${unitPath}.dimension`, str),
        factor: nullable(declaration.factor, `${unitPath}.factor`, num),
        offset: nullable(declaration.offset, `${unitPath}.offset`, num),
      };
    }),
  );
  const unitSets = nullable(raw.unit_sets, `${path}.unit_sets`, (entry, at) =>
    list(entry, at, (set, setPath) => {
      const declaration = obj(set, setPath);
      return {
        id: str(declaration.id, `${setPath}.id`),
        name: str(declaration.name, `${setPath}.name`),
        units: mapOf(declaration.units, `${setPath}.units`, str),
      };
    }),
  );
  const components = nullable(raw.components, `${path}.components`, (entry, at) =>
    list(entry, at, (component, componentPath) => {
      const declaration = obj(component, componentPath);
      return {
        name: str(declaration.name, `${componentPath}.name`),
        molar_mass: nullable(declaration.molar_mass, `${componentPath}.molar_mass`, num),
      };
    }),
  );
  const tools = nullable(raw.tools, `${path}.tools`, (entry, at) =>
    list(entry, at, (tool, toolPath) => {
      const declaration = obj(tool, toolPath);
      return {
        name: str(declaration.name, `${toolPath}.name`),
        description: str(declaration.description, `${toolPath}.description`),
        // A tool's JSON Schema is the agent's contract, carried and not read here.
        input_schema: obj(declaration.input_schema, `${toolPath}.input_schema`),
      };
    }),
  );
  const catalogue: Catalogue = {
    unit_ops: list(raw.unit_ops, `${path}.unit_ops`, form),
  };
  if (units !== null) {
    catalogue.units = units;
  }
  if (unitSets !== null) {
    catalogue.unit_sets = unitSets;
  }
  if (components !== null) {
    catalogue.components = components;
  }
  if (tools !== null) {
    catalogue.tools = tools;
  }
  return catalogue;
}

/** One named field of an object, where the field is all that is wanted and not the shape. */
export function fieldOf(value: unknown, name: string): unknown {
  if (typeof value !== "object" || value === null) {
    return undefined;
  }
  return (value as Record<string, unknown>)[name];
}

/** One of the module's answers, which crosses as a string. */
function parse<T>(json: string, decode: (value: unknown) => T): T {
  let value: unknown;
  try {
    value = JSON.parse(json);
  } catch {
    throw new WireError("the module's answer", "a JSON document", json.slice(0, 60));
  }
  return decode(value);
}

export const parseEnvelope = (json: string): Envelope => parse(json, decodeEnvelope);
export const parseCatalogue = (json: string): Catalogue => parse(json, decodeCatalogue);
