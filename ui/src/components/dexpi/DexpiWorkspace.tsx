/**
 * The Drawing workspace: a DEXPI file in, its coverage and a simulation scaffold out.
 *
 * **A viewer, not a CAD.** It imports a drawing, shows it, and says what azoth would model and
 * what it would not. It does not move a symbol, author one, or edit the file; the drawing is
 * handed to CAD as XML, and the simulation is edited in the Flowsheet workspace.
 */

import { useRef, useState } from "react";

import { saveText } from "../../save";
import type { ExchangeDoor } from "../../wire/session";
import type {
  CoverageStatus,
  DexpiDiagnostic,
  DexpiExport,
  DexpiView,
  ExchangeOptions,
} from "../../wire/types";
import { WorkspaceSwitch, type Workspace } from "../calculations/WorkspaceSwitch";
import { CoveragePanel } from "./CoveragePanel";
import { DrawingSvg } from "./DrawingSvg";

interface Flags {
  allowLossy: boolean;
  projectRecycles: boolean;
  topologyOnly: boolean;
  allowUnconverged: boolean;
}

export function DexpiWorkspace({
  exchange,
  document,
  openable,
  hidden,
  onOpenFlowsheet,
  onWorkspace,
}: {
  exchange: ExchangeDoor | null;
  document: string;
  openable: boolean;
  hidden: boolean;
  onOpenFlowsheet: (toml: string) => boolean;
  onWorkspace: (workspace: Workspace) => void;
}) {
  const [sourceName, setSourceName] = useState<string | null>(null);
  // ui-ok: the DEXPI view is a second document, not part of the flowsheet envelope.
  const [view, setView] = useState<DexpiView | null>(null);
  const [unreadable, setUnreadable] = useState<DexpiDiagnostic[]>([]);
  const [selected, setSelected] = useState<string | null>(null);
  const [filter, setFilter] = useState<CoverageStatus | "all">("all");
  const [showNotDrawn, setShowNotDrawn] = useState(true);
  const [fault, setFault] = useState<string | null>(null);
  const [flags, setFlags] = useState<Flags>({
    allowLossy: false,
    projectRecycles: true,
    topologyOnly: false,
    allowUnconverged: false,
  });
  const [exporting, setExporting] = useState(false);
  const [exportReport, setExportReport] = useState<DexpiExport | null>(null);
  const fileInput = useRef<HTMLInputElement>(null);
  const canvas = useRef<HTMLDivElement>(null);

  const load = async (file: File) => {
    if (exchange === null) {
      setFault("the exchange is not available on this door");
      return;
    }
    setFault(null);
    try {
      const xml = await file.text();
      const response = await exchange.inspect(xml);
      setSourceName(file.name);
      setSelected(null);
      setFilter("all");
      setExportReport(null);
      if (response.ok && response.view !== undefined) {
        setView(response.view);
        setUnreadable([]);
      } else {
        setView(null);
        setUnreadable(response.diagnostics ?? []);
      }
    } catch (error) {
      setFault(String(error));
    }
  };

  const exportProfile = async (profile: "process" | "plant") => {
    if (exchange === null) {
      return;
    }
    const options: ExchangeOptions = {
      allow_lossy: flags.allowLossy,
      project_recycles: profile === "process" ? flags.projectRecycles : false,
      azoth_extensions: true,
      topology_only: profile === "process" ? flags.topologyOnly : false,
      allow_unconverged: profile === "process" ? flags.allowUnconverged : false,
      drawing: true,
      overrides: {},
    };
    setExporting(true);
    try {
      const report = await exchange.exportFlowsheet(document, profile, options);
      setExportReport(report);
      if (report.ok && report.xml !== null) {
        const name = profile === "process" ? "pfd" : "pid";
        await saveText(`${name}.dexpi.xml`, report.xml, "application/xml");
      }
    } catch (error) {
      setFault(String(error));
    } finally {
      setExporting(false);
    }
  };

  const downloadSvg = () => {
    const svg = canvas.current?.querySelector("svg");
    if (svg === null || svg === undefined) {
      return;
    }
    void saveText(`${sourceName ?? "drawing"}.svg`, svg.outerHTML, "image/svg+xml");
  };

  const openFlowsheet = () => {
    if (view === null || view.scaffold === null) {
      return;
    }
    if (onOpenFlowsheet(view.scaffold.toml)) {
      onWorkspace("flowsheet");
    }
  };

  return (
    <div
      className="app dexpi-workspace"
      hidden={hidden}
      onDragOver={(event) => event.preventDefault()}
      onDrop={(event) => {
        event.preventDefault();
        const file = event.dataTransfer.files[0];
        if (file !== undefined) {
          void load(file);
        }
      }}
    >
      <header className="dexpi-head">
        <div className="dexpi-title">
          <span className="dexpi-mark" aria-hidden="true">
            ⌗
          </span>
          <div>
            <strong>{sourceName ?? "Drawing"}</strong>
            <span className="dexpi-sub">
              {view === null
                ? "import a DEXPI 2.0 or Proteus drawing"
                : `${view.format === "proteus" ? "Proteus" : "DEXPI 2.0"} · ${view.model_name}`}
            </span>
          </div>
        </div>
        <WorkspaceSwitch workspace="drawing" onWorkspace={onWorkspace} />
        <div className="dexpi-actions">
          <input
            ref={fileInput}
            type="file"
            accept=".xml,.dexpi,application/xml,text/xml"
            hidden
            onChange={(event) => {
              const file = event.target.files?.[0];
              if (file !== null && file !== undefined) {
                void load(file);
              }
              event.target.value = "";
            }}
          />
          <button type="button" onClick={() => fileInput.current?.click()}>
            Import…
          </button>
          <button
            type="button"
            disabled={!openable || view?.scaffold === null || view === null}
            onClick={openFlowsheet}
          >
            Open as flowsheet
          </button>
          <button
            type="button"
            disabled={view === null || view.scaffold === null || exporting}
            onClick={() => void exportProfile("process")}
          >
            Export PFD
          </button>
          <button
            type="button"
            disabled={view === null || view.scaffold === null || exporting}
            onClick={() => void exportProfile("plant")}
          >
            Export P&ID
          </button>
          <button type="button" disabled={view === null} onClick={downloadSvg}>
            Download SVG
          </button>
        </div>
      </header>

      <div className="dexpi-note">
        Read-only. This is not a CAD editor: hand the XML to CAD for symbols and revisions, and
        edit the simulation in the Flowsheet workspace.
      </div>

      {fault === null ? null : <p className="dexpi-fault">{fault}</p>}

      <div className="dexpi-body">
        <div className="dexpi-figure" ref={canvas}>
          {view?.drawing === null || view === null ? (
            <div className="dexpi-empty-figure">
              <p>Import a DEXPI P&amp;ID or PFD to see the drawing and what azoth would model.</p>
              <button type="button" onClick={() => fileInput.current?.click()}>
                Choose a file
              </button>
              {unreadable.length === 0 ? null : (
                <ul className="dexpi-unreadable">
                  {unreadable.map((diagnostic, index) => (
                    <li key={`${diagnostic.code}-${index}`}>
                      <strong>{diagnostic.code}</strong> {diagnostic.message}
                    </li>
                  ))}
                </ul>
              )}
            </div>
          ) : (
            <DrawingSvg
              drawing={view.drawing}
              view={view}
              selected={selected}
              onSelect={setSelected}
            />
          )}
        </div>

        {view === null ? null : (
          <CoveragePanel
            view={view}
            selected={selected}
            onSelect={setSelected}
            filter={filter}
            onFilter={setFilter}
            showNotDrawn={showNotDrawn}
            onShowNotDrawn={setShowNotDrawn}
          />
        )}
      </div>

      <footer className="dexpi-foot">
        <div className="dexpi-flags">
          <label>
            <input
              type="checkbox"
              checked={flags.allowLossy}
              onChange={(event) => setFlags({ ...flags, allowLossy: event.target.checked })}
            />
            allow lossy mappings
          </label>
          <label>
            <input
              type="checkbox"
              checked={flags.projectRecycles}
              onChange={(event) => setFlags({ ...flags, projectRecycles: event.target.checked })}
            />
            project recycles (PFD)
          </label>
          <label>
            <input
              type="checkbox"
              checked={flags.topologyOnly}
              onChange={(event) => setFlags({ ...flags, topologyOnly: event.target.checked })}
            />
            topology only (PFD)
          </label>
          <label>
            <input
              type="checkbox"
              checked={flags.allowUnconverged}
              onChange={(event) => setFlags({ ...flags, allowUnconverged: event.target.checked })}
            />
            publish unconverged (PFD)
          </label>
        </div>
        {exportReport === null ? null : exportReport.ok ? (
          <p className="dexpi-export-ok">
            Export written: {exportReport.diagram_groups} drawing group(s),{" "}
            {exportReport.diagram_primitives} primitive(s).
          </p>
        ) : (
          <div className="dexpi-export-refused">
            <strong>The export was refused:</strong>
            <ul>
              {exportReport.diagnostics.map((diagnostic, index) => (
                <li key={`${diagnostic.code}-${index}`}>
                  <strong>{diagnostic.code}</strong> {diagnostic.message}
                </li>
              ))}
            </ul>
          </div>
        )}
      </footer>
    </div>
  );
}
