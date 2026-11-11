/**
 * The Cyberpunk HUD's rail switches and the command palette's pointer entry.
 *
 * **They live in the telemetry strip, not in the deck.** A toggle that hid the deck from inside
 * the deck would take its own control away; the strip stays. The labels are short because the bar
 * also carries the case, the unit set and the Appearance control; the full name is the accessible
 * label and the tooltip.
 */

export interface HudControl {
  id: string;
  label: string;
  short: string;
}

export function HudControls({
  controls,
  state,
  onToggle,
  onCommands,
}: {
  controls: readonly HudControl[];
  state: Record<string, boolean>;
  onToggle: (id: string) => void;
  onCommands: () => void;
}) {
  return (
    <div className="hud-controls" role="group" aria-label="HUD panels">
      {controls.map((control) => {
        const shown = state[control.id] !== false;
        return (
          <button
            key={control.id}
            type="button"
            className="hud-toggle"
            aria-pressed={shown}
            aria-label={`${shown ? "Hide" : "Show"} ${control.label}`}
            title={`${shown ? "Hide" : "Show"} ${control.label}`}
            onClick={() => onToggle(control.id)}
          >
            {control.short}
          </button>
        );
      })}
      <button
        type="button"
        className="hud-toggle"
        aria-haspopup="dialog"
        aria-label="Open the command palette"
        title="Open the command palette (Ctrl+K)"
        onClick={onCommands}
      >
        Cmd
      </button>
    </div>
  );
}
