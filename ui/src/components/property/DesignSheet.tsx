import type { Units } from "../../state/units";
import { ParameterField } from "../ParameterField";
import type { EditorCommand } from "../../wire/commands";
import type { Catalogue, Form, GraphNode } from "../../wire/types";

/**
 * A unit operation's declaration, filled in.
 *
 * `middleware.md`'s claim is that a unit-op window *is* its port declaration, and this is the
 * concrete form of it - the fields are the entry's parameters, with the units and the kinds the
 * library declared and the bounds the model states. Nothing here is a second schema: a parameter
 * the entry does not declare cannot be reached, because only declared ones are drawn.
 *
 * **This is the first sheet of an instance's window**, so it is what a person sees the moment they
 * select a node - which is why the wiring it used to end with is now the Connections sheet rather
 * than a heading below the form.
 */
export function DesignSheet({
  catalogue,
  node,
  units,
  onCommand,
}: {
  catalogue: Catalogue | null;
  node: GraphNode;
  /**
   * The unit a reader wants values in, or `null` where there is no catalogue yet.
   *
   * **A field with a unit reads and writes in this set**, which is what makes the switcher a set
   * of units rather than a way of reading: the document keeps its spec's unit and the conversion
   * happens at the field's two edges.
   */
  units: Units | null;
  onCommand: (command: EditorCommand) => void;
}) {
  const form: Form | undefined = catalogue?.unit_ops.find((entry) => entry.id === node.data.unit);
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
    <div className="sheet" id="sheet-design" role="tabpanel" aria-labelledby="tab-design">
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
    </div>
  );
}
