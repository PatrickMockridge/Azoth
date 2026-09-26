/**
 * Which tab the property view is on, and the rule that keeps it legal.
 *
 * **The remembered tab is per *kind* and not per object**, which is HYSYS's behaviour and the
 * reason this is a map rather than one value: select a second column and you land on the tab you
 * were last on for a unit operation, not on the first one over again.
 *
 * **The failure that makes the rule necessary**: a remembered tab the *new* object does not offer -
 * you were reading a column's stage table and select a pump, which has no stages. `tabFor` is where
 * that is decided, once, rather than in a component that would have to notice.
 */

import type { GraphEdge, GraphNode } from "../wire/types";

/** What kind of object the property view is showing. */
export type Kind = "instance" | "feed" | "product" | "edge";

/** One sheet of the property view. */
export type TabId =
  | "design"
  | "conditions"
  | "composition"
  | "connections"
  | "worksheet"
  | "convergence"
  | "stages"
  | "profiles"
  | "performance"
  | "results";

/** What a selection is, in the property view's terms, or `null` when nothing is selected. */
export function kindOf(node: GraphNode | null, edge: GraphEdge | null): Kind | null {
  if (node !== null) {
    return node.role === "instance" ? "instance" : node.role === "input" ? "feed" : "product";
  }
  return edge === null ? null : "edge";
}

/**
 * The sheets a kind offers, in the order they are shown.
 *
 * **Connections is last**, as it is in the unit-op window this is imitating: the design is what a
 * person came for and the wiring is what they check afterwards.
 */
export function tabsFor(
  kind: Kind,
  edge: GraphEdge | null,
  /**
   * The sheets a *result* adds, inserted after the first one.
   *
   * **The order is this file's and the membership is not.** Which result sheets exist depends on
   * what the run published (`state/results.ts` decides), and where they sit is a rule about the
   * window - so a column's stage table appears where a person looks for the operation's own answer,
   * between its design and its wiring, without this module knowing what a stage is.
   */
  extra: readonly TabId[] = [],
): readonly TabId[] {
  switch (kind) {
    case "instance":
      return ["design", ...extra, "worksheet", "connections"];
    case "feed":
    case "product":
      return ["conditions", "composition", "connections"];
    case "edge":
      // **A convergence is a tear's, and a connection is not torn.** The seven settings live on
      // `[[recycles]]`, so offering the sheet for a plain connection would be a tab that could
      // only ever say "there is nothing here".
      return edge?.data.kind === "recycle"
        ? ["connections", "convergence"]
        : ["connections"];
  }
}

/** The tab to draw: the one remembered for this kind, where it is still there to draw. */
export function tabFor(
  kind: Kind,
  remembered: Readonly<Partial<Record<Kind, TabId>>>,
  available: readonly TabId[],
): TabId {
  const wanted = remembered[kind];
  if (wanted !== undefined && available.includes(wanted)) {
    return wanted;
  }
  // Every kind has at least one sheet, so this is a real tab rather than a fallback value.
  return available[0] ?? "design";
}
