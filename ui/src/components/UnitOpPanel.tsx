/**
 * A unit operation's window: its declaration, filled in.
 *
 * `middleware.md`'s claim is that a unit-op window *is* its port declaration, and this is the
 * concrete form of it — the fields are the entry's parameters, with the units and the kinds the
 * library declared and the bounds the model states. Nothing here is a second schema: a parameter
 * the entry does not declare cannot be reached, because only declared ones are drawn.
 */

import { ParameterField } from "./ParameterField";
import type { EditorCommand } from "../wire/commands";
import type { Catalogue, Envelope, Form, GraphNode } from "../wire/types";

export interface UnitOpPanelProps {
  catalogue: Catalogue | null;
  envelope: Envelope;
  node: GraphNode;
  onCommand: (command: EditorCommand) => void;
}

export function UnitOpPanel({ catalogue, envelope, node, onCommand }: UnitOpPanelProps) {
  const form: Form | undefined = catalogue?.unit_ops.find((entry) => entry.id === node.data.unit);
  const values = node.data.parameters ?? {};

  return (
    <>
      <h2>{node.data.unit_name ?? node.data.name}</h2>
      {form === undefined ? (
        <p className="note">
          `{node.data.unit}` is not a unit operation this palette carries, so there is no
          declaration to fill in. The checker says the same thing as an `unknown_unit_op`.
        </p>
      ) : (
        <>
          {form.runnable ? null : (
            <p className="note">
              {form.refusal ?? "the executor refuses this entry"}
            </p>
          )}
          {form.parameters.length === 0 ? (
            <p className="note">this entry declares no parameters</p>
          ) : null}
          {form.parameters.map((parameter) => (
            <ParameterField
              key={parameter.name}
              parameter={parameter}
              value={values[parameter.name]}
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
              the model also bounds {form.unmodelled_ranges.join(", ")}, which no field here
              carries — the run can refuse those and this form cannot show them
            </p>
          )}
          <h2>Where it goes</h2>
          <p className="note">
            {envelope.flowsheet.graph.edges
              .filter((edge) => edge.source === node.id || edge.target === node.id)
              .map((edge) => edge.data.path)
              .join(", ") || "nothing is wired to it yet"}
          </p>
        </>
      )}
    </>
  );
}
