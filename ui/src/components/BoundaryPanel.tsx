/**
 * The boundary: what the flowsheet is fed and what it produces.
 *
 * **A row selects; a separate control removes.** The panel used to make a feed's own row the
 * remove action, so clicking one to look at it deleted it and every connection it fed. Selecting
 * is the common gesture and removing is the rare one, so they are two controls now.
 *
 * **A feed is edited as a record, not as a comma-separated string.** `FeedEditor` owns the
 * component picker, the fraction basis and the one `add_input`; this panel only lists the
 * boundary and decides which feed is being edited.
 */

import { useState } from "react";

import { FeedEditor, newFeedDraft } from "./feed/FeedEditor";
import type { Units } from "../state/units";
import type { Catalogue, Envelope } from "../wire/types";
import type { EditorCommand } from "../wire/commands";

export interface BoundaryPanelProps {
  catalogue: Catalogue | null;
  envelope: Envelope;
  /** The unit a reader wants values in. A feed is stored in the units its spec declares. */
  units: Units;
  onCommand: (command: EditorCommand) => void;
  /** A feed or product row selects the object, as the canvas does. */
  onSelect: (id: string) => void;
}

export function BoundaryPanel({
  catalogue,
  envelope,
  units,
  onCommand,
  onSelect,
}: BoundaryPanelProps) {
  const nodes = envelope.flowsheet.graph.nodes;
  const feeds = nodes.filter((node) => node.role === "input");
  const products = nodes.filter((node) => node.role === "product");

  return (
    <div className="group">
      <div className="label">Boundary</div>

      {feeds.map((feed) => (
        <div className="boundary-row" key={feed.id}>
          <button type="button" className="entry" onClick={() => onSelect(feed.id)}>
            <span className="name">{feed.data.name}</span>
            <span className="tag">feed</span>
          </button>
          <button
            type="button"
            className="remove"
            aria-label={`remove feed ${feed.data.name}`}
            title="remove this feed and every connection that draws from it"
            onClick={() => onCommand({ command: "remove_input", name: feed.data.name })}
          >
            ×
          </button>
        </div>
      ))}
      {products.map((product) => (
        <div className="boundary-row" key={product.id}>
          <button type="button" className="entry" onClick={() => onSelect(product.id)}>
            <span className="name">{product.data.name}</span>
            <span className="tag">product</span>
          </button>
          <button
            type="button"
            className="remove"
            aria-label={`remove product ${product.data.name}`}
            title="remove this product and every connection that writes to it"
            onClick={() => onCommand({ command: "remove_product", name: product.data.name })}
          >
            ×
          </button>
        </div>
      ))}

      <NewFeed
        catalogue={catalogue}
        envelope={envelope}
        units={units}
        onCommand={onCommand}
      />
      <NewProduct envelope={envelope} onCommand={onCommand} />
    </div>
  );
}

/** The record of a feed: a fluid, a flow, a pressure and a temperature. */
function NewFeed({
  catalogue,
  envelope,
  units,
  onCommand,
}: {
  catalogue: Catalogue | null;
  envelope: Envelope;
  units: Units;
  onCommand: (command: EditorCommand) => void;
}) {
  const [open, setOpen] = useState(false);
  const taken = envelope.flowsheet.graph.nodes.map((node) => node.data.name);

  if (!open) {
    return (
      <button type="button" className="entry" onClick={() => setOpen(true)}>
        <span className="name">+ New feed</span>
      </button>
    );
  }

  return (
    <FeedEditor
      catalogue={catalogue}
      units={units}
      initial={newFeedDraft(taken, units)}
      editing={false}
      onApply={(record) => {
        onCommand({ command: "add_input", ...record });
        setOpen(false);
      }}
      onCancel={() => setOpen(false)}
    />
  );
}

/** A product is a name and nothing else: an output is calculated. */
function NewProduct({
  envelope,
  onCommand,
}: {
  envelope: Envelope;
  onCommand: (command: EditorCommand) => void;
}) {
  const suggested = nextName(envelope, "product");
  const [name, setName] = useState(suggested);
  return (
    <div className="field">
      <div className="head">
        <span className="name">+ New product</span>
      </div>
      <input
        value={name}
        aria-label="product name"
        onChange={(event) => setName(event.target.value)}
        placeholder={suggested}
      />
      <button
        type="button"
        className="entry"
        disabled={name.trim() === ""}
        onClick={() => {
          onCommand({ command: "add_product", name: name.trim() });
          setName("");
        }}
      >
        <span className="name">Add the product</span>
      </button>
    </div>
  );
}

/** A name nothing in the graph holds, so a second add does not replace the first. */
function nextName(envelope: Envelope, stem: string): string {
  const taken = new Set(envelope.flowsheet.graph.nodes.map((node) => node.data.name));
  for (let index = 1; index < 100; index += 1) {
    const candidate = `${stem}_${index}`;
    if (!taken.has(candidate)) {
      return candidate;
    }
  }
  return `${stem}_x`;
}
