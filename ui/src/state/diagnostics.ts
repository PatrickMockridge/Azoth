/**
 * What the checker found, grouped by the object it is about and said in the reader's terms.
 *
 * **One row per object, not one per fact.** Adding a heater produced two errors, one for each
 * port, and neither named the heater. The facts are still here; the row is the thing a person
 * acts on, and its heading is the object's own name.
 *
 * **The code and the dotted path are a technical detail.** They are what a developer switches
 * on and what a bug report should carry; they are not the sentence a process engineer reads.
 */

import type { Diagnostic, Envelope, GraphEdge, GraphNode } from "../wire/types";

/** One fact inside a group. */
export interface DiagnosticMessage {
  code: string;
  text: string;
  /** `underfed_port · instances/p1.ports.outlet`, for the disclosure. */
  detail: string;
}

/** One object, and what is wrong with it. */
export interface DiagnosticGroup {
  key: string;
  title: string;
  severity: "error" | "warning";
  /** The node or edge id a click selects, or null where the fact is about no selectable object. */
  select: string | null;
  messages: DiagnosticMessage[];
  /** `not connected yet`, `needs attention` or `warning`; a short word for the heading. */
  summary: string;
}

/** The node a node-like diagnostic points at. */
function nodeFor(target: Diagnostic["target"], envelope: Envelope): GraphNode | undefined {
  switch (target.kind) {
    case "node":
      return envelope.flowsheet.graph.nodes.find(
        (node) => node.role === target.role && node.data.name === target.id,
      );
    case "handle":
    case "parameter": {
      const name = target.kind === "handle" ? target.node : target.node;
      return envelope.flowsheet.graph.nodes.find(
        (node) => node.role === "instance" && node.data.name === name,
      );
    }
    case "endpoint":
      return envelope.flowsheet.graph.nodes.find(
        (node) => node.data.name === target.endpoint,
      );
    default:
      return undefined;
  }
}

/** The name a reader knows the object by. */
function titleOf(node: GraphNode | undefined, fallback: string): string {
  if (node === undefined) {
    return fallback;
  }
  if (node.role === "input") {
    return `Feed ${node.data.name}`;
  }
  if (node.role === "product") {
    return `Product ${node.data.name}`;
  }
  return `${node.data.unit_name ?? "Unit"} ${node.data.name}`;
}

/** The edge a pair names, if the document draws one. */
function edgeFor(from: string, to: string, envelope: Envelope): GraphEdge | undefined {
  return envelope.flowsheet.graph.edges.find(
    (edge) => edge.data.from === from && edge.data.to === to,
  );
}

/** The key, title and selection a diagnostic belongs to. */
function groupOf(
  diagnostic: Diagnostic,
  envelope: Envelope,
): { key: string; title: string; select: string | null } {
  const target = diagnostic.target;
  switch (target.kind) {
    case "node":
    case "handle":
    case "parameter": {
      const node = nodeFor(target, envelope);
      const name = target.kind === "handle" ? target.node : target.kind === "parameter" ? target.node : target.id;
      return {
        key: `node:${name}`,
        title: titleOf(node, name),
        select: node?.id ?? null,
      };
    }
    case "edge": {
      const edge = edgeFor(target.from, target.to, envelope);
      return {
        key: `edge:${target.from}->${target.to}`,
        title: `${target.from} → ${target.to}`,
        select: edge?.id ?? null,
      };
    }
    case "endpoint": {
      const node = nodeFor(target, envelope);
      return {
        key: `endpoint:${target.endpoint}`,
        title: titleOf(node, target.endpoint),
        select: node?.id ?? null,
      };
    }
    case "palette":
      return { key: `palette:${target.id}`, title: `Palette · ${target.id}`, select: null };
    case "document":
      return { key: "document", title: "The flowsheet", select: null };
  }
}

/** One fact as the sentence a reader gets. */
function messageOf(diagnostic: Diagnostic): string {
  const target = diagnostic.target;
  if (diagnostic.code === "underfed_port" || diagnostic.code === "overfed_port") {
    const incoming = diagnostic.detail["incoming"] === true;
    const port = target.kind === "handle" ? target.port : "port";
    const count = Number(diagnostic.detail["count"] ?? 0);
    if (diagnostic.code === "underfed_port") {
      return incoming
        ? `The inlet \`${port}\` is not connected.`
        : `The outlet \`${port}\` is not connected to anything.`;
    }
    return incoming
      ? `The inlet \`${port}\` has ${String(count)} connections; it takes one.`
      : `The outlet \`${port}\` has ${String(count)} connections; it produces one.`;
  }
  if (diagnostic.code === "missing_parameter") {
    const parameter = target.kind === "parameter" ? target.name : "a parameter";
    return `Set \`${parameter}\`; it has no value.`;
  }
  return diagnostic.message;
}

/** The object groups an envelope's diagnostics. */
export function groupDiagnostics(envelope: Envelope): DiagnosticGroup[] {
  const groups = new Map<string, DiagnosticGroup>();
  for (const diagnostic of envelope.diagnostics) {
    const { key, title, select } = groupOf(diagnostic, envelope);
    const severity = diagnostic.severity === "warning" ? "warning" : "error";
    let group = groups.get(key);
    if (group === undefined) {
      group = {
        key,
        title,
        severity,
        select,
        messages: [],
        summary: severity === "warning" ? "warning" : "needs attention",
      };
      groups.set(key, group);
    }
    if (severity === "error") {
      group.severity = "error";
      group.summary = "needs attention";
    }
    group.messages.push({
      code: diagnostic.code,
      text: messageOf(diagnostic),
      detail: `${diagnostic.code} · ${diagnostic.section}/${diagnostic.path}`,
    });
  }

  // **Construction first, then errors, then warnings.** A half-built unit is what a reader is
  // doing; a warning is what waits until they are finished.
  for (const group of groups.values()) {
    if (group.messages.every((message) => message.code === "underfed_port" || message.code === "missing_parameter")) {
      group.summary = "not connected yet";
    }
  }
  return [...groups.values()].sort((left, right) => {
    const rank = (group: DiagnosticGroup) =>
      group.summary === "not connected yet" ? 0 : group.severity === "error" ? 1 : 2;
    return rank(left) - rank(right) || left.title.localeCompare(right.title);
  });
}
