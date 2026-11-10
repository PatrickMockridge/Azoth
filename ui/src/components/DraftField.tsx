/**
 * A value that is sent when the typing stops.
 *
 * **The last control that sent a command per keystroke.** A numeric field wrote on every change,
 * so typing `1e-5` sent `NaN` at the `e` and the command refused it; a blank field sent an unset
 * on the way to a value. The draft is form state: Apply is blur or Enter, Escape throws it away,
 * and a value that is not a number is not sent at all.
 */

import { useState } from "react";

export interface DraftFieldProps {
  /** The committed value as text; empty means unset. */
  value: string;
  kind: "number" | "text";
  placeholder?: string;
  label?: string;
  onCommit: (value: string | number) => void;
}

export function DraftField({ value, kind, placeholder, label, onCommit }: DraftFieldProps) {
  const [draft, setDraft] = useState<string | null>(null);
  const shown = draft ?? value;

  const commit = () => {
    if (draft === null) {
      return;
    }
    if (draft.trim() === "") {
      onCommit("");
      setDraft(null);
      return;
    }
    if (kind === "number") {
      const number = Number(draft);
      if (!Number.isFinite(number)) {
        // An unreadable number is not a command; put the committed value back.
        setDraft(null);
        return;
      }
      onCommit(number);
    } else {
      onCommit(draft);
    }
    setDraft(null);
  };

  return (
    <input
      value={shown}
      inputMode={kind === "number" ? "decimal" : "text"}
      placeholder={placeholder}
      aria-label={label}
      onChange={(event) => setDraft(event.target.value)}
      onBlur={commit}
      onKeyDown={(event) => {
        if (event.key === "Enter") {
          commit();
        } else if (event.key === "Escape") {
          setDraft(null);
        }
      }}
    />
  );
}
