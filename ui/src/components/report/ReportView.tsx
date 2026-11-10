/**
 * The printable case: what was run, what it reached, and what it rests on.
 *
 * **A view over the report model, and a print dialog over the view.** The report is assembled
 * once from the envelope; this draws the tables and the drawing and calls `window.print()`, which
 * is the browser's PDF path and needs no library. Nothing here reaches the engine.
 */

import { useRef } from "react";

import { Pfd } from "../../report/pfd";
import { buildCaseReport } from "../../report/model";
import { saveText } from "../../save";
import type { Units } from "../../state/units";
import type { Catalogue, Envelope } from "../../wire/types";

export interface ReportViewProps {
  envelope: Envelope;
  catalogue: Catalogue | null;
  units: Units;
  onClose: () => void;
}

export function ReportView({ envelope, catalogue, units, onClose }: ReportViewProps) {
  const report = buildCaseReport(envelope, catalogue, units, new Date().toLocaleString());
  const svg = useRef<SVGSVGElement>(null);

  const downloadSvg = () => {
    const drawing = svg.current?.outerHTML;
    if (drawing === undefined) {
      return;
    }
    const stem = report.title.replace(/[^A-Za-z0-9._-]+/g, "_") || "case";
    void saveText(`${stem}-pfd.svg`, drawing, "image/svg+xml");
  };

  return (
    <section className="report" role="dialog" aria-label="Case report">
      <header className="report-head">
        <div>
          <h2>{report.title}</h2>
          <p className="note">
            {report.id} · {report.generated} · {report.status} · units {report.unitSet} ·{" "}
            {report.library} {report.version}
          </p>
        </div>
        <div className="report-actions">
          <button type="button" className="primary" onClick={() => window.print()}>
            Print / PDF
          </button>
          <button type="button" onClick={downloadSvg}>
            Download drawing (SVG)
          </button>
          <button type="button" onClick={onClose}>
            Close
          </button>
        </div>
      </header>

      <p className="report-screen-note note">
        The port is incomplete and every shipped model is marked unverified or partially verified.
        These are screening values, not design deliverables.
      </p>

      <div className="report-pfd">
        <Pfd envelope={envelope} svgRef={svg} />
      </div>

      <p className="report-summary">
        {report.summary.unitOps} unit operation(s) · {report.summary.streams} stream(s) ·{" "}
        {report.summary.warnings} caveat(s) · {report.summary.errors} error(s)
      </p>

      {report.tables.map((table) => (
        <section className="report-table" key={table.title}>
          <h3>{table.title}</h3>
          {table.rows.length === 0 ? (
            <p className="note">nothing to show</p>
          ) : (
            <table className="grid">
              <thead>
                <tr>
                  {table.headers.map((header) => (
                    <th scope="col" key={header}>
                      {header}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {table.rows.map((row, index) => (
                  <tr key={`${table.title}-${String(index)}`}>
                    {row.map((cell, at) => (
                      <td key={`${String(at)}-${cell}`}>{cell}</td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          )}
        </section>
      ))}

      <section className="report-table">
        <h3>{report.diagnostics.title}</h3>
        {report.diagnostics.rows.length === 0 ? (
          <p className="note">the model check has no objection to this document</p>
        ) : (
          <table className="grid">
            <thead>
              <tr>
                {report.diagnostics.headers.map((header) => (
                  <th scope="col" key={header}>
                    {header}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {report.diagnostics.rows.map((row, index) => (
                <tr key={`diagnostic-${String(index)}`}>
                  {row.map((cell, at) => (
                    <td key={`${String(at)}-${cell}`}>{cell}</td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </section>

      <p className="report-statement">
        The report is a reading of one run of the document above. It is not a design deliverable;
        each model and its verification are printed in Model checks.
      </p>
    </section>
  );
}
