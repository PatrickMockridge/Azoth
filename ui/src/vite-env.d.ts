/// <reference types="vite/client" />

/**
 * The commit the static deploy was built from.
 *
 * Injected by `.github/workflows/pages.yml` and absent in a development build, where the title
 * strip shows the library version alone. It is a presentation fact and is not on the wire.
 */
interface ImportMetaEnv {
  readonly VITE_AZOTH_BUILD?: string;
}

/**
 * The shipped flowsheets, inlined at build time.
 *
 * **Read from the repository rather than copied into `public/`.** A copy is a second document that
 * can go stale against the one the tests run, and this import makes the dependency explicit: the
 * editor opens the same `specs/flowsheets/` the CLI runs and the Rust tests hold to the NeqSim
 * capture — the demo it starts on, and the blank one New opens.
 */
declare module "../../specs/flowsheets/demo.toml?raw" {
  const document: string;
  export default document;
}

declare module "../../specs/flowsheets/blank.toml?raw" {
  const document: string;
  export default document;
}
