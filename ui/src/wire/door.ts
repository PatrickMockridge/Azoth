/**
 * Which door the page asked for.
 *
 * **A URL parameter and not a control in the tool bar.** A backend is a property of where the page
 * was loaded from rather than something to toggle while editing, and naming it in the URL is what
 * makes the hosted editor reproducible — as a link someone can be sent, and as a case a test can
 * put the app in without a widget to click.
 *
 * `?serve=http://127.0.0.1:4000` is the hosted door; its absence is the module. The server must be
 * told the editor's own origin (`azoth serve --allow-origin http://localhost:5173`), because a
 * JSON `POST` is not a "simple" request and the answer is unreadable without the CORS header.
 */

import { moduleDoors } from "./client";
import { HttpDoor, HttpExchange, HttpSheetDoor } from "./http";
import type { Doors } from "./session";

/**
 * The served base a query string names, or `null` for the module.
 *
 * A trailing slash is trimmed because every call is `${base}/rpc`; an empty or blank value is the
 * module, because `?serve=` with nothing after it names no server rather than a server at "".
 */
export function serveUrl(search: string): string | null {
  const named = new URLSearchParams(search).get("serve");
  if (named === null) {
    return null;
  }
  const base = named.trim().replace(/\/+$/, "");
  return base === "" ? null : base;
}

/** Both doors this page was loaded for: the flowsheet and the calculation book. */
export async function door(search: string): Promise<Doors> {
  const base = serveUrl(search);
  return base === null
    ? moduleDoors()
    : {
        flowsheet: new HttpDoor(base),
        sheets: new HttpSheetDoor(base),
        exchange: new HttpExchange(base),
      };
}
