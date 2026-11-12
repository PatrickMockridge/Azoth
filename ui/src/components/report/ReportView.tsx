/**
 * The printable case: what was run, what it reached, and what it rests on.
 *
 * **A view over the report model, and a print dialog over the view.** The report is assembled
 * once from the envelope; this draws the tables and the drawing and calls `window.print()`, which
 * is the browser's PDF path and needs no library. Nothing here reaches the engine.
 *
 * **The cover reads, the appendix reproduces.** The results on the cover are plain labels, the
 * reader's units and the library's sentence about each number. The appendix keeps the model's
 * field names, library units and provenance, and the same run can be downloaded as JSON so an
 * agent can reproduce and compare it. One standing notice, not one per table.
 */

import { useEffect, useRef } from "react";

import { Pfd } from "../../report/pfd";
import { buildCaseReport, type ReportTable } from "../../report/model";
import { buildRunRecord } from "../../report/runRecord";
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
  const dialog = useRef<HTMLElement>(null);

  // **A dialog moves focus in and gives it back.** The overlay covered the editor, but the
  // keyboard stayed behind it; Escape and Tab were the editor's.
  useEffect(() => {
    const previous = document.activeElement;
    dialog.current?.querySelector<HTMLElement>("button")?.focus();
    return () => {
      if (previous instanceof HTMLElement) {
        previous.focus();
      }
    };
  }, []);

  useEffect(() => {
    const onKey = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        onClose();
      }
    };
    window.addEventListener("keydown", onKey);
    return () => window.removeEventListener("keydown", onKey);
  }, [onClose]);

  const stem = report.title.replace(/[^A-Za-z0-9._-]+/g, "_") || "case";

  const downloadSvg = () => {
    const drawing = svg.current?.outerHTML;
    if (drawing === undefined) {
      return;
    }
    void saveText(`${stem}-pfd.svg`, drawing, "image/svg+xml");
  };

  const downloadRecord = () => {
    const record = buildRunRecord(envelope, units, new Date().toISOString());
    void saveText(`${stem}-run.json`, JSON.stringify(record, null, 2), "application/json");
  };

  return (
    <section
      ref={dialog}
      className="report"
      role="dialog"
      aria-modal="true"
      aria-label="Case report"
      onKeyDown={(event) => {
        if (event.key !== "Tab") {
          return;
        }
        const buttons = [...(dialog.current?.querySelectorAll<HTMLButtonElement>("button") ?? [])];
        const first = buttons[0];
        const last = buttons[buttons.length - 1];
        if (first === undefined || last === undefined) {
          return;
        }
        if (event.shiftKey && document.activeElement === first) {
          last.focus();
          event.preventDefault();
        } else if (!event.shiftKey && document.activeElement === last) {
          first.focus();
          event.preventDefault();
        }
      }}
    >
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
          <button type="button" onClick={downloadRecord}>
            Download run record (JSON)
          </button>
          <button type="button" onClick={onClose}>
            Close
          </button>
        </div>
      </header>

      <div className="report-pfd">
        <Pfd envelope={envelope} svgRef={svg} />
      </div>

      <p className="report-summary">
        {report.summary.unitOps} unit operation(s) · {report.summary.streams} stream(s) ·{" "}
        {report.summary.warnings} caveat(s) · {report.summary.errors} error(s)
      </p>

      {report.cover.map((table) => (
        <Table key={table.title} table={table} />
      ))}

      <section className="report-appendix">
        <h2>Appendix — machine-readable record</h2>
        <p className="note">
          The model&apos;s own field names and units, the parameters the run was given, and the
          provenance. Download the run record for the same values in JSON.
        </p>
        {report.appendix.map((table) => (
          <Table key={table.title} table={table} />
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
      </section>

      <p className="report-statement">
        Screening values. Every number is a result object with its provenance; the appendix and the
        run record carry the model, the verification status and the hashes for an independent
        cross-check.
      </p>
    </section>
  );
}

/** One table, with its note. An empty table is not printed: a heading over nothing is noise. */
function Table({ table }: { table: ReportTable }) {
  if (table.rows.length === 0) {
    return null;
  }
  return (
    <section className="report-table">
      <h3>{table.title}</h3>
      {table.note === undefined ? null : <p className="note">{table.note}</p>}
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
    </section>
  );
}
