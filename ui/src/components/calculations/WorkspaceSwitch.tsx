/** The two documents a reader can move between; the flowsheet stays loaded either way. */
export type Workspace = "flowsheet" | "calculations";

export function WorkspaceSwitch({
  workspace,
  onWorkspace,
}: {
  workspace: Workspace;
  onWorkspace: (workspace: Workspace) => void;
}) {
  return (
    <div className="workspace-switch" role="group" aria-label="Workspace">
      <button
        type="button"
        className={workspace === "flowsheet" ? "active" : ""}
        aria-pressed={workspace === "flowsheet"}
        onClick={() => onWorkspace("flowsheet")}
      >
        Flowsheet
      </button>
      <button
        type="button"
        className={workspace === "calculations" ? "active" : ""}
        aria-pressed={workspace === "calculations"}
        onClick={() => onWorkspace("calculations")}
      >
        Calculations
      </button>
    </div>
  );
}
