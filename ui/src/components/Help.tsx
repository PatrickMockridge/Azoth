/**
 * The guide, in the dock.
 *
 * **Explanations of this editor, not of the models.** The model's own assumptions live in the
 * model's page in the book; repeating them here would be a second copy that could disagree.
 * Everything below is a convention a newcomer meets in the first ten minutes and would otherwise
 * have to infer from a tooltip.
 */

export function Help() {
  return (
    <div className="help">
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
        <dd>Click a unit operation, a feed, a product or a stream. The navigator on the left does the same.</dd>
        <dt>Move</dt>
        <dd>Drag a unit operation. Nothing is recomputed by moving it.</dd>
        <dt>Connect</dt>
        <dd>Drag from an outlet handle to an inlet handle. A stream is named by the outlet that produces it.</dd>
        <dt>Remove</dt>
        <dd>Select and press Backspace. A removed feed takes its connections with it.</dd>
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
        <strong>Messages</strong> holds the model check's findings and any run error; click one to
        select the object it is about. A unit operation&apos;s <strong>Results</strong> sheet holds
        what it reached that is not on an outlet stream, and its <strong>Caveats</strong> are the
        calculation&apos;s own warnings. <strong>Checks</strong> names the model and how far it is
        verified.
      </p>

      <h2>What this is not</h2>
      <p>
        The port is incomplete and every shipped model is marked unverified or partially
        verified. Treat the numbers as screening values, not as design deliverables; the book&apos;s
        front page states what is missing.
      </p>
    </div>
  );
}
