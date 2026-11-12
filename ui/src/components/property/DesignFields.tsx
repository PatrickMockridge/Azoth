/**
 * A unit operation's declared parameters, as controls.
 *
 * **The Inspector and the node window are two surfaces over one writer.** They both render this
 * component, so the field mapping belongs here once: the Inspector stacks the fields because its
 * column is narrow, and the node window uses the row layout because it has the width to put the
 * library's explanation beside the control. Neither is a second editor.
 */

import { ParameterField } from "../ParameterField";
import type { Units } from "../../state/units";
import type { EditorCommand } from "../../wire/commands";
import type { Catalogue, GraphNode } from "../../wire/types";

export function DesignFields({
  catalogue,
  node,
  units,
  onCommand,
  layout = "stacked",
}: {
  catalogue: Catalogue | null;
  node: GraphNode;
  /** The unit a reader wants values in, or `null` where there is no catalogue yet. */
  units: Units | null;
  onCommand: (command: EditorCommand) => void;
  /** Stacked is the Inspector's narrow column; row is the node window's wider sheet. */
  layout?: "stacked" | "row";
}) {
  const form = catalogue?.unit_ops.find((entry) => entry.id === node.data.unit);
  const values = node.data.parameters ?? {};

  if (form === undefined) {
    return (
      <p className="note">
        `{node.data.unit}` is not a unit operation this palette carries, so there is no declaration
        to fill in. The checker says the same thing as an `unknown_unit_op`.
      </p>
    );
  }

  return (
    <>
      {form.runnable ? null : <p className="note">{form.refusal ?? "the executor refuses this entry"}</p>}
      {form.parameters.length === 0 ? (
        <p className="note">this entry declares no parameters</p>
      ) : null}
      {form.parameters.map((parameter) => (
        <ParameterField
          key={parameter.name}
          parameter={parameter}
          value={values[parameter.name]}
          units={units}
          layout={layout}
          onChange={(value) =>
            onCommand(
              value === undefined || value === null || value === ""
                ? {
                    command: "unset_parameter",
                    instance: node.data.name,
                    name: parameter.name,
                  }
                : {
                    command: "set_parameter",
                    instance: node.data.name,
                    name: parameter.name,
                    value,
                  },
            )
          }
        />
      ))}
      {form.unmodelled_ranges.length === 0 ? null : (
        <p className="note">
          the model also bounds {form.unmodelled_ranges.join(", ")}, which no field here carries —
          the run can refuse those and this form cannot show them
        </p>
      )}
    </>
  );
}
