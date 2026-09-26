import { useState } from "react";

import { kindOf, tabFor, tabsFor, type Kind, type TabId } from "../../state/layout";
import type { EditorCommand } from "../../wire/commands";
import type { Catalogue, Envelope, Form, GraphEdge, GraphNode } from "../../wire/types";
import { CompositionSheet } from "./CompositionSheet";
import { ConditionsSheet } from "./ConditionsSheet";
import { ConnectionsSheet } from "./ConnectionsSheet";
import { DesignSheet } from "./DesignSheet";
import { DockTabs } from "./DockTabs";
import { TearSheet } from "./TearSheet";
import { WorksheetSheet } from "./WorksheetSheet";

/** What each sheet is called on its tab. */
const LABELS: Record<TabId, string> = {
  design: "Design",
  conditions: "Conditions",
  composition: "Composition",
  connections: "Connections",
  worksheet: "Worksheet",
  convergence: "Convergence",
  stages: "Stages",
  profiles: "Profiles",
  performance: "Performance",
  results: "Results",
};

/**
 * The selected object's window, in the shape a process simulator's is: a header, a row of sheets,
 * and the sheet.
 *
 * **The remembered sheet is per *kind***, so selecting a second unit operation lands on the sheet
 * you were last on for a unit operation rather than on Design again — see `state/layout.ts` for
 * the rule that keeps it legal when the new object does not offer it.
 */
export function PropertyView({
  catalogue,
  envelope,
  node,
  edge,
  onCommand,
  onSelect,
}: {
  catalogue: Catalogue | null;
  envelope: Envelope;
  node: GraphNode | null;
  edge: GraphEdge | null;
  onCommand: (command: EditorCommand) => void;
  onSelect: (id: string | null) => void;
}) {
  const [remembered, setRemembered] = useState<Partial<Record<Kind, TabId>>>({});
  const kind = kindOf(node, edge);

  if (kind === null) {
    return (
      <p className="note">
        nothing is selected — click a unit operation, a stream or a connection
      </p>
    );
  }

  const available = tabsFor(kind, edge);
  const active = tabFor(kind, remembered, available);
  const form: Form | undefined =
    node === null ? undefined : catalogue?.unit_ops.find((entry) => entry.id === node.data.unit);

  return (
    <>
      <div className="prop-header">
        <h2>{title(kind, node, edge)}</h2>
        <p className="note">{subtitle(kind, node, edge, form)}</p>
      </div>
      <DockTabs
        label={title(kind, node, edge)}
        active={active}
        tabs={available.map((id) => ({ id, label: LABELS[id] }))}
        onTab={(id) => setRemembered({ ...remembered, [kind]: id as TabId })}
      />
      {sheet(active, { catalogue, envelope, node, edge, onCommand, onSelect })}
    </>
  );
}

/** The window's heading: what this object is called. */
function title(kind: Kind, node: GraphNode | null, edge: GraphEdge | null): string {
  switch (kind) {
    case "instance":
      return node?.data.unit_name ?? node?.data.name ?? "unit operation";
    case "feed":
      return "Feed";
    case "product":
      return "Product";
    case "edge":
      return edge?.data.kind === "recycle" ? `Tear ${edge.data.path}` : "Connection";
  }
}

/** The line under it: which object it is, and anything about it the heading cannot say. */
function subtitle(
  kind: Kind,
  node: GraphNode | null,
  edge: GraphEdge | null,
  form: Form | undefined,
): string {
  if (kind === "edge") {
    return `${edge?.data.from ?? ""} → ${edge?.data.to ?? ""}`;
  }
  const name = node?.data.name ?? "";
  if (kind === "instance" && form !== undefined && !form.runnable) {
    return `${name} · the executor refuses this entry`;
  }
  if (kind === "instance" && form === undefined) {
    return `${name} · not a unit operation this palette carries`;
  }
  return name;
}

/** The sheet itself, by id. */
function sheet(
  active: TabId,
  context: {
    catalogue: Catalogue | null;
    envelope: Envelope;
    node: GraphNode | null;
    edge: GraphEdge | null;
    onCommand: (command: EditorCommand) => void;
    onSelect: (id: string | null) => void;
  },
) {
  const { catalogue, envelope, node, edge, onCommand, onSelect } = context;
  switch (active) {
    case "design":
      return node === null ? null : (
        <DesignSheet catalogue={catalogue} node={node} onCommand={onCommand} />
      );
    case "conditions":
      return node === null ? null : (
        <ConditionsSheet envelope={envelope} node={node} onCommand={onCommand} />
      );
    case "composition":
      return node === null ? null : <CompositionSheet envelope={envelope} node={node} />;
    case "worksheet":
      return node === null ? null : <WorksheetSheet envelope={envelope} node={node} />;
    case "convergence":
      return edge === null ? null : (
        <TearSheet catalogue={catalogue} edge={edge} onCommand={onCommand} />
      );
    case "connections":
      return (
        <ConnectionsSheet
          envelope={envelope}
          node={node}
          edge={edge}
          onCommand={onCommand}
          onSelect={onSelect}
        />
      );
    default:
      // The result sheets arrive with the workbook, which is where the results region is read.
      return <p className="note">there is nothing on this sheet yet</p>;
  }
}
