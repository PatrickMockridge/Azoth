/**
 * A result field's name as a person reads it.
 *
 * **Presentation only, and never the record.** The model's own field name is the canonical one:
 * it is what the case file, the Python record and the JSON run record carry, and it is what the
 * machine-readable appendix prints. This module turns that name into a label for the human
 * summary - `outlet_t` is a temperature, not an "outlet t" - using a small token map for the
 * state variables and the acronyms this library uses, and snake_case-to-words for everything
 * else. It is deliberately not a per-field table: a field added to a model next year still gets
 * a readable label, and no second naming of the library's own fields exists here.
 */

/** The tokens whose bare spelling is not the word a reader needs. */
const TOKENS: Record<string, string> = {
  n: "molar flow",
  p: "pressure",
  t: "temperature",
  h: "molar enthalpy",
  z: "composition",
};

/** Acronyms that stay upper-case rather than becoming a word. */
const ACRONYMS = new Set([
  "io",
  "lmtd",
  "npsh",
  "nps",
  "sil",
  "thd",
  "ua",
  "vfd",
  "fs",
  "hv",
  "lv",
  "co2",
  "h2s",
]);

/**
 * A field name as a label: `outlet_duty` → "Outlet duty", `outlet_t` → "Outlet temperature".
 *
 * The first word is capitalised and the rest are left as the model wrote them, except that an
 * acronym is upper-cased and a state variable becomes its quantity's name.
 */
export function humanLabel(key: string): string {
  const words = key
    .split("_")
    .filter((word) => word !== "")
    .map((word) => {
      const lower = word.toLowerCase();
      const token = TOKENS[lower];
      if (token !== undefined) {
        return token;
      }
      return ACRONYMS.has(lower) ? lower.toUpperCase() : word;
    });
  if (words.length === 0) {
    return key;
  }
  const [first, ...rest] = words;
  return [first === undefined ? key : capitalise(first), ...rest].join(" ");
}

function capitalise(word: string): string {
  return word.charAt(0).toUpperCase() + word.slice(1);
}

/**
 * The library's sentence as a plain explanation.
 *
 * The spec writes field names in backticks - `` `inlet_n * (outlet_h - inlet_h)` `` - which is
 * right for a technical document and noise beside a value in the editor and the report. The
 * formula stays; only the quoting goes.
 */
export function explanationOf(description: string | null | undefined): string | null {
  if (description === null || description === undefined) {
    return null;
  }
  const text = description.replaceAll("`", "").trim();
  return text === "" ? null : text;
}
