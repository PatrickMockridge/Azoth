/**
 * Search the databank instead of spelling a substance.
 *
 * **The catalogue is the authority on names, so this offers names and nothing else.** A free-text
 * entry would let a typo travel to Solve and come back as a run error, which is a worse place to
 * learn that the fluid's name is `n-butane` and not `nbutane`.
 */

import { useState } from "react";

import type { Catalogue } from "../../wire/types";

/** How many matches the picker shows at once; the search is meant to narrow, not to browse. */
const LIMIT = 8;

export interface ComponentPickerProps {
  catalogue: Catalogue | null;
  /** Names already in the feed; a duplicate is not offered. */
  chosen: string[];
  onPick: (name: string) => void;
}

export function ComponentPicker({ catalogue, chosen, onPick }: ComponentPickerProps) {
  const [query, setQuery] = useState("");
  const needle = query.trim().toLowerCase();
  const taken = new Set(chosen);
  const matches =
    needle === ""
      ? []
      : (catalogue?.components ?? [])
          .filter(
            (component) =>
              !taken.has(component.name) && component.name.toLowerCase().includes(needle),
          )
          .slice(0, LIMIT);

  return (
    <div className="component-picker">
      <input
        value={query}
        onChange={(event) => setQuery(event.target.value)}
        placeholder="Search substances, e.g. methane"
        aria-label="Search substances"
      />
      {needle === "" ? (
        <p className="hint">Type a few letters of the substance&apos;s name.</p>
      ) : matches.length === 0 ? (
        <p className="hint bad">no substance in the databank matches that</p>
      ) : (
        <div className="matches">
          {matches.map((component) => (
            <button
              key={component.name}
              type="button"
              onClick={() => {
                onPick(component.name);
                setQuery("");
              }}
            >
              {component.name}
            </button>
          ))}
        </div>
      )}
    </div>
  );
}
