/**
 * What a click selects, from a node id or from a diagnostic's target.
 *
 * **A diagnostic carries a target precisely so this file is a reading and not a guess.** The
 * checker derives one from the variant's own fields, so "put a mark on this" is the library's
 * answer; the only work here is turning it into the single node a canvas can select.
 */

import type { Envelope, GraphNode, Target } from "../wire/types";

/** The node a selection names, or none. */
export function selectedNode(envelope: Envelope, id: string | null): GraphNode | null {
  if (id === null) {
    return null;
  }
  return envelope.flowsheet.graph.nodes.find((node) => node.id === id) ?? null;
}

/**
 * The node to select for a diagnostic, where one can be said.
 *
 * `null` is a real answer: an endpoint that names no instance, a palette entry and a fact about
 * the document whole are none of them about a node, and selecting something arbitrary to look
 * responsive would be worse than selecting nothing.
 */
export function targetNodeId(target: Target, envelope: Envelope): string | null {
  switch (target.kind) {
    case "node":
      return `${target.role}:${target.id}`;
    case "handle":
      return `${target.role}:${target.node}`;
    case "parameter":
      return `instance:${target.node}`;
    case "edge": {
      // An edge is about two endpoints; the producing one is the node a user would go to.
      const producer = target.from.split(".")[0] ?? "";
      const id = `instance:${producer}`;
      return envelope.flowsheet.graph.nodes.some((node) => node.id === id) ? id : null;
    }
    default:
      return null;
  }
}
