/// <reference types="vite/client" />

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
