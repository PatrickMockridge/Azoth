/**
 * What a blank page needs, in the order a flowsheet is built.
 *
 * Shown only when the document has no objects. It is a hint rather than a wizard: every step is
 * a control that already exists on the screen, and the text says where it is.
 */
export function StarterHints() {
  return (
    <div className="starter-hints">
      <h2>An empty flowsheet</h2>
      <ol>
        <li>
          Add a <strong>feed</strong> from the Boundary panel on the left: a fluid, a flow, a
          pressure and a temperature.
        </li>
        <li>
          Click a unit operation in the <strong>palette</strong> to place it, then fill in its
          Design sheet.
        </li>
        <li>
          Drag from the feed&apos;s handle to the unit&apos;s inlet, and from its outlet to a{" "}
          <strong>product</strong>.
        </li>
        <li>
          Press <strong>Solve</strong>, and read the streams in the Workbook.
        </li>
      </ol>
      <p className="note">
        The <strong>Case</strong> menu can open a worked example instead, and the example is a
        faster way to see what a finished document looks like.
      </p>
    </div>
  );
}
