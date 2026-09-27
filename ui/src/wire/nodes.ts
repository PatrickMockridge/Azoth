/**
 * A node id, as the command that removes it.
 *
 * **The TypeScript reading of `split_node_id`.** A projection's node id is `{role}:{name}`, and the
 * role decides which command removes the thing: an instance goes by `remove_instance`, a feed and
 * a product by their own. This is a wire mapping and not a selection, which is why it lives here
 * and not in `state/selection.ts`.
 *
 * It returns `null` for an id whose role is none of the three rather than guessing one: a node the
 * projection did not make is not a node to invent a command for, and a right-click that silently
 * did nothing to the wrong thing is worse than one that does nothing.
 */

import type { EditorCommand } from "./commands";
import { asRole } from "./decode";

/** The command that removes the node an id names, or none. */
export function removeCommandFor(nodeId: string): EditorCommand | null {
  const separator = nodeId.indexOf(":");
  if (separator < 0) {
    return null;
  }
  // **Read, not asserted.** The id's prefix is a projection convention, so a role is one of three
  // strings or it is none of them — and `as Role` would have made a fourth look like a first.
  const role = asRole(nodeId.slice(0, separator));
  if (role === null) {
    return null;
  }
  const name = nodeId.slice(separator + 1);
  // **Three cases and no `default`.** The type has three inhabitants once `asRole` has read the
  // prefix, so a fourth branch could not be reached and a `default: return null` would be a case
  // that says the reader is not exhaustive when it is. A role added to the projection fails
  // `tsc` here rather than silently removing nothing.
  switch (role) {
    case "instance":
      return { command: "remove_instance", id: name };
    case "input":
      return { command: "remove_input", name };
    case "product":
      return { command: "remove_product", name };
  }
}
