/// <reference types="vite/client" />

/**
 * The shipped flowsheet, inlined at build time.
 *
 * **Read from the repository rather than copied into `public/`.** A copy is a second document that
 * can go stale against the one the tests run, and this import makes the dependency explicit: the
 * editor opens the same `specs/flowsheets/demo.toml` the CLI runs and the Rust tests hold to the
 * NeqSim capture.
 */
declare module "../../specs/flowsheets/demo.toml?raw" {
  const document: string;
  export default document;
}
