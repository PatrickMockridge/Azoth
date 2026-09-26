/**
 * The flowsheet's objects, as a list of things to select.
 *
 * **A navigator lists the document, not the library.** Every row is an object this flowsheet
 * actually has — a unit op, a stream, a feed, a product — so a document with one node has one row
 * and a blank one has none. The palette beside it is the opposite: it lists what *can* be added,
 * which no document changes.
 *
 * **A row's id is the id the canvas selects by** — `instance:hx1` for a node, `e5` for an edge —
 * because two panes that select one thing must be one selection rather than two that agree. Nothing
 * here is a second selection model, and nothing here reads a value: the panel is a way somewhere.
 *
 * **Families are ordered as the catalogue declares them**, so this pane and the palette put
 * `two port` in the same place, and rows within a family keep the document's own order. That is the
 * one thing this module decides rather than reads.
 */

import type { Catalogue, Envelope, Form, GraphNode } from "../wire/types";

/** One row: a thing to select, what it is called, and what to say about it. */
export interface NavRow {
  /** The id the canvas selects by, which is the whole point of the row. */
  id: string;
  label: string;
  /** A word to the right of the label: `recycle`, or nothing. */
  tag: string | null;
}

/** One heading and the rows under it. */
export interface NavGroup {
  label: string;
  rows: NavRow[];
}

/**
 * The family a form was filed under, or `unfiled` where the catalogue did not say.
 *
 * **`unfiled` is a visible group rather than a silent merge into `other`**: an entry that arrives
 * without its family is a wire that stopped carrying it, and a person should see that as its own
 * heading rather than as one more entry in a group that looks like a family. The palette groups by
 * the same word, which is why it is one function rather than two strings.
 */
export function familyOfForm(form: Form): string {
  return form.family ?? "unfiled";
}

/**
 * A family as a heading rather than as a spec directory: `two_port` is two words to a reader.
 *
 * Beside `familyOfForm` because the palette and this pane must spell it the same way — one is the
 * name a group is keyed by and the other is what it says, and the pair drifts if it lives twice.
 */
export function familyLabel(family: string): string {
  return family.replace("_", " ");
}

/** The family a unit op instance belongs to, read through the catalogue's entry for it. */
function familyOfNode(node: GraphNode, catalogue: Catalogue | null): string {
  const form = (catalogue?.unit_ops ?? []).find((entry) => entry.id === node.data.unit);
  return form === undefined ? "unfiled" : familyOfForm(form);
}

/**
 * The whole tree: unit operations by family, then the streams, then the boundary.
 *
 * The three boundary groups are last because they are the ends of the flowsheet rather than its
 * middle, and a recycle appears among the streams with a word beside it — a tear is a stream the
 * solver iterates, not a separate kind of object.
 */
export function navigationOf(envelope: Envelope, catalogue: Catalogue | null): NavGroup[] {
  const { nodes, edges } = envelope.flowsheet.graph;

  const byFamily = new Map<string, NavRow[]>();
  const feeds: NavRow[] = [];
  const products: NavRow[] = [];

  for (const node of nodes) {
    if (node.role === "instance") {
      const family = familyOfNode(node, catalogue);
      byFamily.set(family, [
        ...(byFamily.get(family) ?? []),
        { id: node.id, label: node.data.name, tag: null },
      ]);
    } else if (node.role === "input") {
      feeds.push({ id: node.id, label: node.data.name, tag: null });
    } else {
      products.push({ id: node.id, label: node.data.name, tag: null });
    }
  }

  const streams = edges.map((edge) => ({
    id: edge.id,
    label: edge.data.path,
    tag: edge.data.kind === "recycle" ? "recycle" : null,
  }));

  const groups: NavGroup[] = [];
  for (const family of familyOrder(catalogue, byFamily)) {
    groups.push({ label: familyLabel(family), rows: byFamily.get(family) ?? [] });
  }
  if (streams.length > 0) {
    groups.push({ label: "Streams", rows: streams });
  }
  if (feeds.length > 0) {
    groups.push({ label: "Feeds", rows: feeds });
  }
  if (products.length > 0) {
    groups.push({ label: "Products", rows: products });
  }
  return groups;
}

/**
 * The families to draw, in the catalogue's order, with anything it did not carry after them.
 *
 * A family the catalogue declares but this document does not use is not a group: an empty heading
 * would say the document has one, and the navigator is a reading of the document.
 */
function familyOrder(catalogue: Catalogue | null, byFamily: Map<string, NavRow[]>): string[] {
  const declared = (catalogue?.unit_ops ?? []).map(familyOfForm);
  const ordered: string[] = [];
  for (const name of declared) {
    if (name !== "unfiled" && byFamily.has(name) && !ordered.includes(name)) {
      ordered.push(name);
    }
  }
  for (const name of byFamily.keys()) {
    if (!ordered.includes(name)) {
      ordered.push(name);
    }
  }
  return ordered;
}
