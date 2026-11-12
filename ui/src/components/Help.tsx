/**
 * The guide, in the dock.
 *
 * **Explanations of this editor, not of the models.** The model's own assumptions live in the
 * model's page in the book; repeating them here would be a second copy that could disagree.
 * Everything below is a convention a newcomer meets in the first ten minutes and would otherwise
 * have to infer from a tooltip.
 */

import { manualHome } from "../manual";
import { definitionOf } from "../state/verification";

/** The four statuses the library derives, in the order its own enum lists them. */
const VERIFICATION_STATUSES = [
  "verified",
  "partially_verified",
  "source_needed",
  "unverified",
] as const;

export function Help() {
  return (
    <div className="help">
      <p>
        <a href={manualHome()} target="_blank" rel="noreferrer">
          Open the user manual
        </a>{" "}
        for the task guides, the calculation reference and the full validity notes.
      </p>

      <h2>Getting started</h2>
      <p>
        A flowsheet is <strong>feeds</strong>, <strong>unit operations</strong> and{" "}
        <strong>products</strong>, joined by streams. Add a unit operation from the palette, click
        it to fill in its Design sheet, then connect it by dragging from an outlet handle to an
        inlet handle. Press <strong>Solve</strong> when the document is complete.
      </p>
      <p>
        The status beside the title says what the document is: <strong>refused</strong> means the
        model check has found something to fix, <strong>stale</strong> means it has not been run
        since the last edit, and <strong>solved</strong> means the numbers are the current run's.
      </p>

      <h2>The drawing</h2>
      <dl>
        <dt>Select</dt>
        <dd>Click a unit operation, a feed, a product or a stream. The navigator on the left does the same. Double-click a node for its window: a unit operation&apos;s Design fields are editable there, and below them are what it is wired to and what it reached.</dd>
        <dt>Move</dt>
        <dd>Drag a unit operation. Nothing is recomputed by moving it.</dd>
        <dt>Connect</dt>
        <dd>Drag from an outlet handle to an inlet handle. A stream is named by the outlet that produces it. Every connection point is labelled <strong>IN</strong> or <strong>OUT</strong>; a two-stream exchanger labels <strong>hot</strong> and <strong>cold</strong> on each side. <strong>One line uses one point</strong>: a mixer has a separate inlet for every feed, and a splitter a separate outlet for every product. Connectors are drawn as right-angle routes.</dd>
        <dt>Remove</dt>
        <dd>Select and press Delete or Backspace, or use Remove in the object&apos;s window. A removed feed takes its connections with it. Undo restores it.</dd>
        <dt>Tear</dt>
        <dd>A connection can be turned into a recycle from its Connections sheet. A tear is drawn moving and is iterated until its residual is small.</dd>
      </dl>

      <h2>Feeds and products</h2>
      <p>
        A feed is a fluid and four numbers: molar flow, pressure, temperature and composition.
        Use the component picker; the databank decides whether a substance name exists. A product
        is a name only, because its stream is calculated.
      </p>

      <h2>Units</h2>
      <p>
        The unit set at the top right changes what is shown; it never changes the document. A
        field reads and accepts the set in force, so a temperature typed in °F arrives in the
        model as the kelvin it stands for.
      </p>

      <h2>Messages and results</h2>
      <p>
        <strong>Messages</strong> holds the model check&apos;s findings and any run error; click one
        to select the object it is about. A unit operation&apos;s <strong>Results</strong> sheet
        holds what it reached that is not on an outlet stream: a plain label, the value and unit,
        and the library&apos;s own sentence about how the number is obtained — often the equation
        behind it, so the answer can be checked by hand. Only fields the run published are drawn,
        and its <strong>Caveats</strong> are the calculation&apos;s own warnings. The collapsed
        <strong> Machine-readable fields</strong> block and the <strong>Provenance</strong> tab
        hold the same values in the model&apos;s field names and the hashes behind them.
      </p>

      <h2>Provenance</h2>
      <p>
        The <strong>Provenance</strong> tab answers what a number rests on: the model, the NeqSim
        class it was ported from, the data authority (the standard keycard this build ships), the
        spec and kernel hashes, and the model&apos;s verification status. The status says how far
        the model itself is checked; it is not the status of the run, which <strong>Report</strong>
        and the run record carry separately.
      </p>
      <dl>
        {VERIFICATION_STATUSES.map((status) => (
          <div key={status}>
            <dt>
              <code>{status}</code>
            </dt>
            <dd>{definitionOf(status).definition}</dd>
          </div>
        ))}
      </dl>

      <h2>Calculation sheets</h2>
      <p>
        The <strong>Calculations</strong> switch opens a calculation book beside the flowsheet.
        Pick a calculation from the list, fill its inputs and press <strong>Calculate</strong>;
        <strong>Load worked example</strong> fills the calculation&apos;s own published example.
        A book holds several sheets, is saved in TOML, and can be printed as a report. The four
        calculations whose inputs or outputs are tables are listed and not yet available. See
        the{" "}
        <a
          href={`${manualHome()}using/calculations.html`}
          target="_blank"
          rel="noreferrer"
        >
          calculation-sheet guide
        </a>{" "}
        for the layout and the unit-set behaviour.
      </p>

      <h2>What this is not</h2>
      <p>
        The port is incomplete and every shipped model is marked unverified or partially
        verified. Treat the numbers as screening values, not as design deliverables;{" "}
        <a
          href={`${manualHome()}index.html#not-for-design-work-yet`}
          target="_blank"
          rel="noreferrer"
        >
          the manual&apos;s front page
        </a>{" "}
        states what is missing.
      </p>
    </div>
  );
}
