import { seriesOf } from "../../state/results";
import type { Units } from "../../state/units";
import type { UnitResult } from "../../wire/types";

/** How many points a sparkline draws before it is no longer a glance. */
const MAX_POINTS = 64;

/**
 * The vectors a unit operation filled in: a column's trays, a reactor's axial profile, a Gibbs
 * solver's energy history.
 *
 * **A table first, and a drawing beside it.** The table is the data and the surface a screen reader
 * walks; the sparkline is what makes a divergence visible at a glance, which is the whole reason a
 * process engineer looks at a tray profile. A sheet with only the chart would be a picture of the
 * answer with no way to read it.
 *
 * **The chart is drawn from the same numbers and computes nothing**: the x axis is the index, and
 * each series is scaled to its own range because two profiles of different magnitudes share a
 * sheet.
 */
export function StagesSheet({
  result,
  units,
}: {
  result: UnitResult;
  /** The unit a reader wants the values in; a catalogue that has not loaded converts nothing. */
  units: Units;
}) {
  const series = seriesOf(result, units);

  if (series.length === 0) {
    return (
      <p className="note">
        this unit operation published no profile — a mixer's answer is its stream, and only a stage,
        a reactor or an iterative solve fills a vector
      </p>
    );
  }

  const rows = Math.max(...series.map((entry) => entry.values.length));

  return (
    <div className="sheet" id="sheet-stages" role="tabpanel" aria-labelledby="tab-stages">
      <table className="grid">
        <thead>
          <tr>
            <th scope="col">#</th>
            {series.map((entry) => (
              <th key={entry.key} scope="col" className="num">
                {entry.key}
                {entry.unit === null ? null : <span className="unit"> {entry.unit}</span>}
              </th>
            ))}
          </tr>
        </thead>
        <tbody>
          {Array.from({ length: rows }, (_, index) => (
            <tr key={index}>
              <th scope="row">{index}</th>
              {series.map((entry) => (
                <td key={entry.key} data-num="true">
                  {entry.values[index] === undefined
                    ? "—"
                    : entry.values[index]?.toFixed(5)}
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      {series.map((entry) => (
        <div className="spark" key={entry.key}>
          <div className="spark-label">
            {entry.key}
            {entry.unit === null ? null : <span className="unit"> {entry.unit}</span>}
          </div>
          <Sparkline values={entry.values} />
        </div>
      ))}
    </div>
  );
}

/**
 * One series as a line, in SVG.
 *
 * **`vector-effect="non-scaling-stroke"`**, because a scaled `viewBox` multiplies the stroke width
 * and a one-pixel profile line becomes a thick one on a wide panel - the reason a hand-rolled chart
 * is subtly wrong if that is left out.
 *
 * A `null` would break the line rather than plunge to zero, which is what a missing tray value
 * means; the wire carries no nulls inside a vector today, so the guard is for the day one does.
 */
export function Sparkline({ values }: { values: readonly number[] }) {
  if (values.length < 2) {
    return <p className="note">one point is not a profile</p>;
  }
  const drawn = values.length > MAX_POINTS ? values.slice(0, MAX_POINTS) : values;
  const low = Math.min(...drawn);
  const high = Math.max(...drawn);
  // A flat series has no range to scale into, so it is drawn down the middle rather than at nought.
  const span = high - low === 0 ? 1 : high - low;
  const height = 28;
  const width = 240;
  const points = drawn
    .map((value, index) => {
      const x = (index / (drawn.length - 1)) * width;
      const y = height - ((value - low) / span) * height;
      return `${x.toFixed(2)},${y.toFixed(2)}`;
    })
    .join(" ");

  return (
    <svg
      className="chart"
      viewBox={`0 0 ${String(width)} ${String(height)}`}
      width={width}
      height={height}
      role="img"
      aria-label={`${String(drawn.length)} points from ${low.toPrecision(4)} to ${high.toPrecision(4)}`}
    >
      <polyline points={points} vectorEffect="non-scaling-stroke" />
    </svg>
  );
}
