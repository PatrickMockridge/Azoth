/**
 * The editor's appearance: a layout and a colour, and the two places that fact lives.
 *
 * **The default is CSS's, not this module's.** `:root` carries the Windows XP light tokens, so a
 * document with no attributes and no script already has the default look - which is what makes it
 * a default rather than a value somebody set. This module only writes an override.
 *
 * **Only this module touches the attributes or the store at runtime.** `index.html` also sets them
 * before the first paint so a returning reader never sees the default flash; that script reads the
 * same key and the same four values, and `ui/test/appearance.test.ts` holds the two copies equal.
 */

/** The colour axis: brightness and ink, within a layout. */
export type Theme = "light" | "neutral" | "dark";

/** The layout axis: the chrome's lineage. */
export type Layout = "xp" | "mac" | "cyberpunk";

/** One point in the two-dimensional choice. */
export interface Appearance {
  layout: Layout;
  theme: Theme;
}

/** What a first visit gets: the software the audience most likely came from. */
export const DEFAULT: Appearance = { layout: "xp", theme: "light" };

/** The one key this module writes. */
const KEY = "azoth.appearance";

/** The key a pre-appearance build wrote, read once so a returning reader keeps their choice. */
const LEGACY = "azoth.theme";

function isLayout(value: unknown): value is Layout {
  return value === "xp" || value === "mac" || value === "cyberpunk";
}

function isTheme(value: unknown): value is Theme {
  return value === "light" || value === "neutral" || value === "dark";
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === "object" && value !== null;
}

/** A stored value, or an empty object when the store is absent or holds something else. */
function stored(): Record<string, unknown> {
  try {
    const raw = window.localStorage.getItem(KEY);
    if (raw === null) {
      return {};
    }
    const parsed: unknown = JSON.parse(raw);
    return isRecord(parsed) ? parsed : {};
  } catch {
    // A browser can refuse storage (a private window, a policy), and a corrupt value is not a
    // fault worth a message: the default is still a complete appearance.
    return {};
  }
}

/** The URL's preview choice, where it names a valid layout or colour. */
function preview(): Partial<Appearance> {
  const params = new URLSearchParams(window.location.search);
  const layout = params.get("layout");
  const theme = params.get("theme");
  return {
    ...(isLayout(layout) ? { layout } : {}),
    ...(isTheme(theme) ? { theme } : {}),
  };
}

/** Whether this page load is a preview, which is not remembered unless the reader chooses. */
export function isPreview(): boolean {
  const chosen = preview();
  return chosen.layout !== undefined || chosen.theme !== undefined;
}

/** The appearance in force: the URL's preview, else what was chosen, else the default. */
export function readAppearance(): Appearance {
  const remembered = stored();
  let theme: Theme = DEFAULT.theme;
  if (isTheme(remembered.theme)) {
    theme = remembered.theme;
  } else {
    try {
      // An older build remembered only dark or light. The layout that build had is XP's ancestor.
      const legacy = window.localStorage.getItem(LEGACY);
      if (isTheme(legacy)) {
        theme = legacy;
      }
    } catch {
      // As above: the default stands.
    }
  }
  return {
    layout: isLayout(remembered.layout) ? remembered.layout : DEFAULT.layout,
    theme,
    ...preview(),
  };
}

/**
 * Put an appearance in force, and remember it.
 *
 * `persist: false` is the pre-paint effect of a URL preview: the page shows what the link asked
 * for, but a reader who did not choose it has not changed their saved preference.
 */
export function applyAppearance(
  appearance: Appearance,
  { persist = true }: { persist?: boolean } = {},
): void {
  document.documentElement.dataset.layout = appearance.layout;
  document.documentElement.dataset.theme = appearance.theme;
  if (persist) {
    try {
      window.localStorage.setItem(KEY, JSON.stringify(appearance));
    } catch {
      // The attributes are the appearance; the store is only how it survives a reload.
    }
  }
  const meta = document.querySelector('meta[name="theme-color"]');
  if (meta !== null) {
    const app = getComputedStyle(document.documentElement).getPropertyValue("--bg-app").trim();
    if (app !== "") {
      meta.setAttribute("content", app);
    }
  }
}
