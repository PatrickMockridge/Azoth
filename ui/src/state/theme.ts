/**
 * The theme: which one is on, and the two places that fact lives.
 *
 * **Dark is the default and it is CSS's, not this module's.** `:root` carries the dark tokens, so
 * a document with no `data-theme` and no script is already dark - which is what makes the default
 * a default rather than a value somebody set. This module only ever *overrides* it.
 *
 * **Only here touches the attribute or the store.** A second reader or writer would be a second
 * chance for the attribute and the remembered choice to disagree, and the disagreement shows as a
 * theme that flips back on the next reload.
 */

/** The two themes a person can choose. */
export type Theme = "dark" | "light";

/** Where the choice is remembered. One key, and nothing else in the app writes one. */
const KEY = "azoth.theme";

/** The theme in force: what was chosen, or dark where nothing was. */
export function readTheme(): Theme {
  try {
    return window.localStorage.getItem(KEY) === "light" ? "light" : "dark";
  } catch {
    // A browser can refuse storage (a private window, a policy). The default is still dark, and a
    // theme that cannot be remembered is not a fault worth a message.
    return "dark";
  }
}

/** Put a theme in force, and remember it. */
export function applyTheme(theme: Theme): void {
  document.documentElement.dataset.theme = theme;
  try {
    window.localStorage.setItem(KEY, theme);
  } catch {
    // As above: the attribute is the theme, and the store is only how it survives a reload.
  }
}
