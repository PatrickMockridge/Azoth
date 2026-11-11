/**
 * Whether the first-run card is still owed.
 *
 * **A reader's fact, like the theme and the unit set**, so it is remembered in `localStorage`
 * and travels with no document. A browser that refuses storage shows the card again rather than
 * failing to load; that is the same trade the unit set makes.
 */

const KEY = "azoth.welcome";

/** The card is shown until the reader has seen and dismissed or used it. */
export function readWelcome(): boolean {
  try {
    return window.localStorage.getItem(KEY) === null;
  } catch {
    return true;
  }
}

/** Record that the card has been dealt with, so the next visit lands on the editor. */
export function rememberWelcome(): void {
  try {
    window.localStorage.setItem(KEY, "seen");
  } catch {
    // Storage can be refused; the card is dismissed for this visit regardless.
  }
}
