//! The official DEXPI XML schema, vendored and pinned.
//!
//! The schema is CC BY 4.0 (DEXPI e.V.), not AGPL: it is distributed unmodified, with its
//! source and digest recorded in `NOTICE`. The bundled file declares `version="1.0.0"` in its
//! own header even though it is the 2.0 release, so **the digest is the identity** - the
//! generator refuses a file whose SHA-256 is not [`SHA256`].
//!
//! Validation itself is `xmllint`'s: a Rust XSD validator is a large dependency and libxml2
//! through FFI is forbidden by the workspace's `unsafe_code`. The test gate fails, rather than
//! skips, when `xmllint` is absent.

/// The schema, byte for byte as published for DEXPI XML V2.0.0 (2.0.1 is identical).
pub const SCHEMA: &str = include_str!("../schema/DEXPI_XML_Schema.xsd");

/// SHA-256 of [`SCHEMA`], as published with the DEXPI 2.0 release.
pub const SHA256: &str = "0bd6568b7c52fd59ac4ef6d8547ab2b2aa890c009a3a1bfbc5d7cc18c53b36e6";
