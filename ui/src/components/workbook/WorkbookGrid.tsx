import type { ReactNode } from "react";

/**
 * The one grid: a worksheet, a composition sheet, a stage table and the workbook are all this.
 *
 * **A real `<table>`**, with `<th scope="col">` and `<th scope="row">`, because a grid is the one
 * surface here a screen reader walks properly - `role="grid"` on a pile of divs would be a claim
 * about keyboard navigation this does not implement, and a table's own semantics are free.
 *
 * **A cell whose run is older than the document is marked**, which is what `dirty` means: the
 * numbers are the previous document's, and a table that showed them plainly would be the one place
 * the editor lied about it.
 */
export interface GridColumn {
  key: string;
  label: string;
  /** The unit, where the column carries one - drawn in the header and not beside every cell. */
  unit?: string | null;
  /** Whether the cells are numbers, which right-aligns them. */
  numeric?: boolean;
}

export function WorkbookGrid({
  columns,
  rows,
  stale,
  onSelect,
  empty,
}: {
  columns: readonly GridColumn[];
  /** One row: the key of its header cell and the cells, positionally with `columns`. */
  rows: readonly { id: string; head: ReactNode; cells: ReactNode[]; title?: string }[];
  /** Whether the values are older than the document, which is `envelope.dirty`. */
  stale: boolean;
  onSelect?: (id: string) => void;
  /** What to say where there are no rows at all, rather than an empty table. */
  empty?: ReactNode;
}) {
  if (rows.length === 0) {
    return <p className="note">{empty ?? "nothing to show"}</p>;
  }

  return (
    <table className="grid">
      <thead>
        <tr>
          {columns.map((column) => (
            <th key={column.key} scope="col" className={column.numeric === true ? "num" : ""}>
              {column.label}
              {column.unit === undefined || column.unit === null ? null : (
                <span className="unit"> {column.unit}</span>
              )}
            </th>
          ))}
        </tr>
      </thead>
      <tbody>
        {rows.map((row) => (
          <tr key={row.id} data-id={row.id}>
            <th scope="row">
              {onSelect === undefined ? (
                row.head
              ) : (
                <button type="button" className="link" onClick={() => onSelect(row.id)}>
                  {row.head}
                </button>
              )}
            </th>
            {row.cells.map((cell, index) => (
              <td
                key={columns[index]?.key ?? index}
                data-num={columns[index]?.numeric === true ? "true" : undefined}
                data-stale={stale ? "true" : undefined}
                title={row.title}
              >
                {cell}
              </td>
            ))}
          </tr>
        ))}
      </tbody>
    </table>
  );
}
