/**
 * The boundary: what the flowsheet is fed and what it produces.
 *
 * **This is where `add_input` becomes reachable for a *new* record.** The inputs panel could edit
 * an existing feed's flow, pressure and temperature, and `components` and `z` were read-only with
 * a note saying the command that changes them is `add_input` — which no control offered. A feed is
 * a record, and the record is what this form writes.
 *
 * **Two rules are enforced here, and both are the checker's own.** Add is disabled until the fluid
 * names at least one substance and the mole fractions match it one for one — which is exactly the
 * two facts `Diagnostic::InputRecord` tests, so the form cannot send a record the checker would
 * refuse for a reason it could see. **Substance names are not checked**: whether `methane` exists
 * is the databank's answer, it arrives as the run's own failure, and guessing at it here would be
 * a second opinion about the databank.
 */

import { useState } from "react";

import type { EditorCommand } from "../wire/commands";
import type { Envelope } from "../wire/types";

export interface BoundaryPanelProps {
  envelope: Envelope;
  onCommand: (command: EditorCommand) => void;
}

export function BoundaryPanel({ envelope, onCommand }: BoundaryPanelProps) {
  const nodes = envelope.flowsheet.graph.nodes;
  const feeds = nodes.filter((node) => node.role === "input");
  const products = nodes.filter((node) => node.role === "product");

  return (
    <div className="group">
      <div className="label">Boundary</div>

      {feeds.map((feed) => (
        <button
          key={feed.id}
          type="button"
          className="entry"
          onClick={() => onCommand({ command: "remove_input", name: feed.data.name })}
        >
          <span className="name">{feed.data.name}</span>
          <span className="tag">feed · remove</span>
        </button>
      ))}
      {products.map((product) => (
        <button
          key={product.id}
          type="button"
          className="entry"
          onClick={() => onCommand({ command: "remove_product", name: product.data.name })}
        >
          <span className="name">{product.data.name}</span>
          <span className="tag">product · remove</span>
        </button>
      ))}

      <NewFeed envelope={envelope} onCommand={onCommand} />
      <NewProduct envelope={envelope} onCommand={onCommand} />
    </div>
  );
}

/** The record of a feed, which is four numbers and a fluid. */
function NewFeed({ envelope, onCommand }: BoundaryPanelProps) {
  const existing = envelope.flowsheet.graph.nodes.find((node) => node.role === "input");
  const [open, setOpen] = useState(false);
  const [name, setName] = useState(nextName(envelope, "feed"));
  const [components, setComponents] = useState(
    (existing?.data.input?.components ?? []).join(", "),
  );
  const [z, setZ] = useState((existing?.data.input?.z ?? []).join(", "));
  const [n, setN] = useState(String(existing?.data.input?.n ?? 1));
  const [P, setP] = useState(String(existing?.data.input?.P ?? 5e5));
  const [T, setT] = useState(String(existing?.data.input?.T ?? 300));

  const fluids = components
    .split(",")
    .map((entry) => entry.trim())
    .filter((entry) => entry !== "");
  const fractions = z
    .split(",")
    .map((entry) => entry.trim())
    .filter((entry) => entry !== "")
    .map(Number);
  // The checker's own two facts, so the form cannot send a record it would refuse.
  const ready =
    name.trim() !== "" &&
    fluids.length >= 1 &&
    fractions.length === fluids.length &&
    fractions.every((fraction) => Number.isFinite(fraction));

  if (!open) {
    return (
      <button type="button" className="entry" onClick={() => setOpen(true)}>
        <span className="name">+ New feed</span>
      </button>
    );
  }

  return (
    <div className="field">
      <div className="head">
        <span className="name">New feed</span>
      </div>
      <input value={name} onChange={(event) => setName(event.target.value)} placeholder="name" />
      <input
        value={components}
        onChange={(event) => setComponents(event.target.value)}
        placeholder="methane, n-butane"
      />
      <input
        value={z}
        onChange={(event) => setZ(event.target.value)}
        placeholder="0.9, 0.1"
      />
      <input value={n} onChange={(event) => setN(event.target.value)} placeholder="mol/s" />
      <input value={P} onChange={(event) => setP(event.target.value)} placeholder="Pa" />
      <input value={T} onChange={(event) => setT(event.target.value)} placeholder="K" />
      {ready ? null : (
        <div className="hint bad">
          a fluid names at least one substance, and one mole fraction each
        </div>
      )}
      <button
        type="button"
        className="entry"
        disabled={!ready}
        onClick={() => {
          onCommand({
            command: "add_input",
            name: name.trim(),
            components: fluids,
            n: Number(n),
            z: fractions,
            P: Number(P),
            T: Number(T),
          });
          setOpen(false);
        }}
      >
        <span className="name">Add the feed</span>
      </button>
    </div>
  );
}

/** A product is a name and nothing else: an output is calculated. */
function NewProduct({ envelope, onCommand }: BoundaryPanelProps) {
  const [name, setName] = useState(nextName(envelope, "product"));
  return (
    <div className="field">
      <div className="head">
        <span className="name">+ New product</span>
      </div>
      <input value={name} onChange={(event) => setName(event.target.value)} placeholder="name" />
      <button
        type="button"
        className="entry"
        disabled={name.trim() === ""}
        onClick={() => onCommand({ command: "add_product", name: name.trim() })}
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
