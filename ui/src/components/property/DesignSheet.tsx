import type { Units } from "../../state/units";
import type { EditorCommand } from "../../wire/commands";
import type { Catalogue, GraphNode } from "../../wire/types";
import { DesignFields } from "./DesignFields";

/**
 * A unit operation's declaration, filled in.
 *
 * **The wrapper is the Inspector's**; the fields and their command mapping are `DesignFields`,
 * shared with the node window so the two surfaces cannot disagree about a parameter.
 *
 * `middleware.md`'s claim is that a unit-op window *is* its port declaration, and this is the
 * concrete form of it - the fields are the entry's parameters, with the units and the kinds the
 * library declared and the bounds the model states. Nothing here is a second schema: a parameter
 * the entry does not declare cannot be reached, because only declared ones are drawn.
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
  return (
    <div className="sheet" id="sheet-design" role="tabpanel" aria-labelledby="tab-design">
      <DesignFields
        catalogue={catalogue}
        node={node}
        units={units}
        onCommand={onCommand}
      />
    </div>
  );
}
