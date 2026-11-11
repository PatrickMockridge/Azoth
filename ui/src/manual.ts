/**
 * Where the user manual lives.
 *
 * **The book is published beside the editor**, at `/azoth/manual/` on Pages, so the production
 * link is relative and survives a different host or a preview deployment. In development the
 * book is not served by Vite, so the link points at the published copy; `VITE_AZOTH_MANUAL`
 * overrides both when a preview deployment has its own book.
 *
 * A calculation's page is generated at `<namespace>/<leaf>.html`, the same path
 * `tools/gen_docs.py` writes.
 */

const configured = (import.meta.env.VITE_AZOTH_MANUAL ?? "").trim();
const base = (
  configured !== ""
    ? configured
    : import.meta.env.DEV
      ? "https://patrickmockridge.github.io/azoth/manual"
      : "manual"
).replace(/\/$/, "");

export function manualHome(): string {
  return `${base}/`;
}

export function calcManualUrl(id: string): string {
  return `${base}/${id.replace(".", "/")}.html`;
}
