//! `azoth mcp` — the tool schema, served over stdio.
//!
//! `docs/src/architecture/middleware.md` says an MCP server is "a projection of this schema rather
//! than a second one", and this is where that is either true or false. **It is true by
//! construction**: the tools are `middleware::tools`'s own list with one key renamed, every call
//! is one `Command` applied to one `Workspace`, and every answer is the envelope. There is no
//! second description of what an edit is anywhere in this file.
//!
//! **Two protocol revisions, because the protocol moved and clients have not all followed.** The
//! current revision, `2026-07-28`, removed the `initialize` handshake: the version travels in each
//! request's `_meta`, `server/discover` is a mandatory RPC, and every result carries `resultType`.
//! `2025-11-25` and earlier are the handshake revisions. The specification says clients and
//! servers **may** support several at once, and a server speaking one lane is unreachable from
//! half its callers — so both are served, and **one encoder** serves them: `resultType` and the
//! cache fields are always written, because the current revision requires them and an older client
//! ignores a field it does not know. What differs is which opening a client uses, and the method
//! set each lane has.
//!
//! **A modern request carries its metadata, and a missing one is malformed.** The current
//! revision's schema puts `protocolVersion` and `clientCapabilities` in `_meta` and requires
//! `_meta` on every request; the specification says a request missing a required field is
//! malformed and **must** be answered `-32602`, with HTTP `400` on that transport. So the lane is
//! entered by what a request carries, with one exception: an `initialize` request selects legacy
//! semantics *"scoped to the stdio process"*, and that scope is why [`Protocol`] exists — a
//! `2025-11-25` client's later `tools/list` carries no `_meta`, because that revision put the
//! version in the handshake, and a server judging every request by its metadata alone would refuse
//! a legacy client's own call as malformed.
//!
//! **The transport is newline-delimited JSON with no embedded newlines.** `serde_json::to_string`
//! emits none, and every message is written and flushed as one line — a client blocks forever on a
//! partial one, which is the classic hand-rolled-framing failure and the reason the test reads a
//! line per request rather than writing everything and reading to the end.
//!
//! **`crate::mcp_http` is this over HTTP**, and it is the same `Protocol` and the same answers: the
//! only things that differ are the ones HTTP owns — status codes, required headers, and the
//! `initialize` refusal a modern-only transport owes a legacy client.
//!
//! **Hand-rolled, and the cost is named.** Five things have to be right: the framing above; error
//! codes (`-32700` parse, `-32601` method, `-32602` a malformed request or an unknown tool,
//! `-32022` the protocol version, which is the specification's own number for it); the required
//! `_meta` fields and which lane a request is in; the `input_schema` → `inputSchema` rename; and
//! the *shape* of every message, which is checked rather than trusted - the specification's JSON
//! Schema for each revision is vendored under `python/tests/validation/vendor/mcp/` and every
//! message this file writes is validated against the definition for its own lane. **What that
//! cannot see is a new revision**: the vendored files are the pin, so upstream moving is invisible
//! until somebody re-fetches them deliberately.

use std::io::{BufRead, Write};
use std::path::Path;

use azoth_core::{AzothError, Result};
use azoth_process::middleware::tools;
use serde_json::{Value, json};

use crate::session::Session;

/// The protocol revisions this server speaks, newest first.
///
/// **The newest is what a client that asks is told**, and both are what a version refusal lists.
/// Read from the specification rather than remembered: `2026-07-28` is current and handshake-free,
/// `2025-11-25` is the last handshake revision.
pub const PROTOCOL_VERSIONS: [&str; 2] = ["2026-07-28", "2025-11-25"];

/// JSON-RPC's own codes, plus the specification's one.
const PARSE_ERROR: i64 = -32700;

/// **`pub(crate)` because `crate::mcp_http` maps it to a status**, and only there: the code is the
/// protocol's and the status is the transport's, so the one place they meet is the one place that
/// needs both.
pub(crate) const METHOD_NOT_FOUND: i64 = -32601;

/// JSON-RPC's `Invalid params`, which the specification also uses for a *malformed request*: a
/// modern request missing a required `_meta` field is `-32602` rather than a code of its own.
const INVALID_PARAMS: i64 = -32602;
const UNSUPPORTED_VERSION: i64 = -32022;

/// The two `_meta` keys `2026-07-28` requires on every request.
const META_VERSION: &str = "io.modelcontextprotocol/protocolVersion";
const META_CAPABILITIES: &str = "io.modelcontextprotocol/clientCapabilities";

/// The methods the current revision defines, for a refusal to name.
const MODERN_METHODS: [&str; 3] = ["server/discover", "tools/list", "tools/call"];

/// The methods the handshake revisions define. **`ping` is here and not above**: the schema for
/// `2026-07-28` defines no `PingRequest`, so a modern client has `server/discover` instead.
const LEGACY_METHODS: [&str; 4] = ["initialize", "ping", "tools/list", "tools/call"];

/// Serve the tool schema over stdio, against one document.
///
/// **One session, held across calls.** The server opens the document once and every `tools/call`
/// edits that same `Workspace`, which is what makes an agent's second call see its first one's
/// effect — and what makes "read the document" unnecessary, since the envelope an edit answers
/// with *is* the document.
///
/// **The document on disk is never written**, which is `azoth edit`'s own property for its own
/// reason: a session is not a file editor. It costs nothing, because `flowsheet.document` is in
/// every envelope, so a client that wants the result writes the string it already has — and that
/// is why `save` is not a tool.
///
/// # Errors
/// A palette that does not load, or a document that does not read. Everything after that is a
/// JSON-RPC error *in the stream* rather than a failure of the server: a client that sent one bad
/// call has not ended the session.
pub fn serve(
    flowsheet: &Path,
    palette_dir: &Path,
    run_after_each_call: bool,
    input: impl BufRead,
    mut output: impl Write,
) -> Result<()> {
    // **One session, and this transport does not know what an edit is.** Everything above the
    // framing is `Session`'s; what is left here is JSON-RPC.
    let mut session = Session::open(flowsheet, palette_dir, run_after_each_call)?;
    let mut protocol = Protocol::new();

    for line in input.lines() {
        let line = line.map_err(|error| AzothError::invalid_input("stdin", error.to_string()))?;
        if line.trim().is_empty() {
            continue;
        }
        let response = match serde_json::from_str::<Value>(&line) {
            Ok(request) => protocol.handle(&mut session, &request),
            // A line that is not JSON is the one failure with no id to answer to, which is why the
            // spec gives it a code of its own and `null` for the id.
            Err(error) => Some(error_response(
                &Value::Null,
                PARSE_ERROR,
                &error.to_string(),
                None,
            )),
        };
        if let Some(response) = response {
            let text = serde_json::to_string(&response)
                .map_err(|error| AzothError::invalid_input("response", error.to_string()))?;
            // One line, flushed: a client reading a line at a time blocks on a partial one.
            writeln!(output, "{text}")
                .and_then(|()| output.flush())
                .map_err(|error| AzothError::invalid_input("stdout", error.to_string()))?;
        }
    }
    Ok(())
}

/// Which revision's semantics a request is served under.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Lane {
    /// The `2026-07-28` revision, where every request carries its own protocol metadata.
    Modern,
    /// The handshake revisions. Entered by `initialize`, and only by `initialize`.
    Legacy,
}

/// The protocol state a connection carries, which is the lane it opened in.
///
/// **A dual-era server selects its behaviour from how the client opens**, in the specification's
/// own words: *"A request carrying modern per-request `_meta` is served statelessly according to
/// this revision. An `initialize` request selects legacy semantics, scoped to the stdio process."*
/// This type is that scope, and it is one struct because both transports need the same answers:
/// `azoth mcp` holds one for the life of the process, and `crate::mcp_http` makes one per request,
/// because HTTP is modern-only and nothing may be assumed from a previous request there.
///
/// **It holds no session**, which is what lets `azoth serve` hand the *same* `Workspace` to this
/// and to `/rpc`: one document, two doors, one encoder.
pub struct Protocol {
    lane: Lane,
}

impl Default for Protocol {
    fn default() -> Self {
        Self::new()
    }
}

impl Protocol {
    /// A connection that has not opened yet, which is the modern lane: a request that carries no
    /// `_meta` and is not a handshake is a *modern request missing its metadata*, which the schema
    /// calls malformed and the specification answers `-32602`.
    #[must_use]
    pub fn new() -> Self {
        Self { lane: Lane::Modern }
    }

    /// One request, as a response — or `None` for a notification.
    ///
    /// **A notification has no reply**, and answering one puts a message on the wire the client is
    /// not reading: the specification's own rule, and the reason this returns an `Option`. On HTTP
    /// that `None` is a `202 Accepted` rather than a silence.
    pub fn handle(&mut self, session: &mut Session, request: &Value) -> Option<Value> {
        // A request with no id is a notification, whatever its method: no reply, because a reply
        // would be a message the client is not reading. This `?` is that rule, spelled here as the
        // absence of an answer; `is_notification` is the same rule named, because the HTTP
        // transport has to ask the question *before* it dispatches rather than after.
        let id = request.get("id").cloned()?;

        let method = request
            .get("method")
            .and_then(Value::as_str)
            .unwrap_or_default();

        if method == "initialize" {
            self.lane = Lane::Legacy;
            return Some(self.initialize(&id, request));
        }
        // Carrying `_meta` is what makes a request modern, whatever came before it on this
        // connection — the specification allows a dual-era server to serve both eras concurrently,
        // so a client that opens with the handshake and then sends metadata is served either way.
        if request.pointer("/params/_meta").is_some() {
            self.lane = Lane::Modern;
        }
        Some(match self.lane {
            Lane::Modern => self.modern(session, &id, method, request),
            Lane::Legacy => self.legacy(session, &id, method, request),
        })
    }

    /// The handshake a `2025-11-25` client opens with, answered whatever came before it.
    ///
    /// The version a client asks for is the one it is told, because this server speaks both — and
    /// one it does not speak is refused with the list, which is the specification's own payload for
    /// that refusal.
    fn initialize(&self, id: &Value, request: &Value) -> Value {
        let asked = request
            .pointer("/params/protocolVersion")
            .and_then(Value::as_str)
            .unwrap_or(PROTOCOL_VERSIONS[0]);
        if !PROTOCOL_VERSIONS.contains(&asked) {
            return error_response(
                id,
                UNSUPPORTED_VERSION,
                &format!("`{asked}` is not a protocol revision this server speaks"),
                Some(unsupported(asked)),
            );
        }
        result(
            id,
            json!({
                "protocolVersion": asked,
                "capabilities": { "tools": { "listChanged": false } },
                "serverInfo": server_info(),
                "instructions": INSTRUCTIONS,
            }),
        )
    }

    /// The current revision's lane: every request carries its own metadata, so nothing is
    /// remembered and nothing may be assumed from a request that came before it.
    ///
    /// **The metadata is required before the version is judged**, because a request with no
    /// `_meta` at all has no version to judge — and the specification gives that case its own
    /// answer (`-32602`, a malformed request) rather than the version refusal's.
    fn modern(&self, session: &mut Session, id: &Value, method: &str, request: &Value) -> Value {
        let version = match meta_version(request) {
            Ok(version) => version,
            Err(sentence) => return error_response(id, INVALID_PARAMS, &sentence, None),
        };
        if !PROTOCOL_VERSIONS.contains(&version) {
            return error_response(
                id,
                UNSUPPORTED_VERSION,
                &format!("`{version}` is not a protocol revision this server speaks"),
                Some(unsupported(version)),
            );
        }

        match method {
            // The current revision's mandatory discovery RPC.
            "server/discover" => result(
                id,
                json!({
                    "supportedVersions": PROTOCOL_VERSIONS,
                    "capabilities": { "tools": { "listChanged": false } },
                    "_meta": { "io.modelcontextprotocol/serverInfo": server_info() },
                    "instructions": INSTRUCTIONS,
                    "ttlMs": 3_600_000,
                    // A property of the software rather than of a document.
                    "cacheScope": "public",
                }),
            ),
            "tools/list" => tools_list(session, id),
            "tools/call" => tools_call(session, id, request),
            // `ping` is not one of them, which is the schema's answer rather than a preference:
            // `2026-07-28` defines no `PingRequest`, and `server/discover` is what replaced it.
            other => unknown_method(id, other, &MODERN_METHODS),
        }
    }

    /// The handshake lane: no request carries `_meta`, because the version was negotiated once.
    fn legacy(&self, session: &mut Session, id: &Value, method: &str, request: &Value) -> Value {
        match method {
            "ping" => result(id, json!({})),
            "tools/list" => tools_list(session, id),
            "tools/call" => tools_call(session, id, request),
            other => unknown_method(id, other, &LEGACY_METHODS),
        }
    }
}

/// The version a modern request declares, or the sentence saying which part is missing.
///
/// **Both keys are required**, which is the schema's own `required` array on its request-metadata
/// object rather than a reading of the prose. `clientCapabilities` is required and *unused*: this
/// server needs no client capability to run a flowsheet — no sampling, no elicitation, no roots —
/// so `-32021 MissingRequiredClientCapability` is never emitted, and the field is checked because
/// the schema says a request without it is malformed.
fn meta_version(request: &Value) -> std::result::Result<&str, String> {
    let Some(meta) = request.pointer("/params/_meta").and_then(Value::as_object) else {
        return Err("a request in this revision carries `params._meta`".to_string());
    };
    let Some(version) = meta.get(META_VERSION).and_then(Value::as_str) else {
        return Err(format!("`_meta` names no `{META_VERSION}`"));
    };
    if !meta.contains_key(META_CAPABILITIES) {
        return Err(format!("`_meta` names no `{META_CAPABILITIES}`"));
    }
    Ok(version)
}

/// The specification's own payload for a version refusal.
///
/// **Both keys, and the schema is why**: its `UnsupportedProtocolVersionError` requires `requested`
/// and `supported` on the error's `data`, so a refusal naming only the versions this server speaks
/// leaves a client unable to say which of its requests was refused. The two sites that emit this
/// are the handshake and the per-request metadata, and they are one function because they are one
/// refusal.
fn unsupported(version: &str) -> Value {
    json!({ "requested": version, "supported": PROTOCOL_VERSIONS })
}

/// A method this lane does not serve — and **the list is the lane's own**, which is the one place
/// the two revisions' method sets are written down beside each other.
fn unknown_method(id: &Value, method: &str, methods: &[&str]) -> Value {
    error_response(
        id,
        METHOD_NOT_FOUND,
        &format!("`{method}` is not a method this server serves"),
        Some(json!({ "methods": methods })),
    )
}

/// The catalogue, which is the same list in both lanes.
fn tools_list(session: &Session, id: &Value) -> Value {
    let list: Vec<Value> = session.tools().iter().map(tool_json).collect();
    result(
        id,
        json!({
            "tools": list,
            "ttlMs": 3_600_000,
            // A property of this build and this palette, so a client may cache it for
            // the session but not across sessions.
            "cacheScope": "private",
        }),
    )
}

/// One tool call, in the MCP result shape — **the same in both lanes**, because a tool call is a
/// tool call: what the revision changed is the metadata around it.
///
/// **The call itself is [`Session::call`]'s**, and this only wraps the answer: the envelope
/// as text for a model to read and as structured content for a client to address, or the
/// library's own sentence when the call is what was wrong. The three refusal states the
/// session keeps apart arrive here as they are - an `isError` result for a call that could not
/// take effect, and `ok: false` inside the envelope for a document the checker refuses.
fn tools_call(session: &mut Session, id: &Value, request: &Value) -> Value {
    let Some(name) = request.pointer("/params/name").and_then(Value::as_str) else {
        return error_response(id, INVALID_PARAMS, "a tool call must name a tool", None);
    };
    // **The spec's own shape for this is a protocol error, not a result**, which is why the
    // name is checked here as well as in the session: a name that is not a tool is `-32602` to
    // an MCP client and a `400` to an HTTP one, and the *list* is one list
    // (`Session::tools`) while the expression differs by transport.
    if !session.tools().iter().any(|tool| tool.name == name) {
        return error_response(id, INVALID_PARAMS, &format!("Unknown tool: {name}"), None);
    }
    let arguments = request
        .pointer("/params/arguments")
        .cloned()
        .unwrap_or_else(|| json!({}));

    match session.call(name, &arguments, None) {
        Ok(envelope) => {
            let text = session
                .envelope_json()
                .unwrap_or_else(|_| envelope.to_string());
            result(
                id,
                json!({
                    "content": [ { "type": "text", "text": text } ],
                    // The same document twice: once for a model to read as text and once for a
                    // client to address field by field. One encoder wrote both.
                    "structuredContent": envelope,
                    "isError": false,
                }),
            )
        }
        Err(sentence) => refused(id, &sentence),
    }
}

/// A call that could not take effect, as an `isError` result rather than a protocol error.
///
/// The specification puts input-validation failures here rather than in a JSON-RPC error, so
/// that a model can read the sentence and correct itself — and the sentence is the library's
/// own, which names the offending field.
fn refused(id: &Value, reason: &str) -> Value {
    result(
        id,
        json!({
            "content": [ { "type": "text", "text": reason } ],
            "isError": true,
        }),
    )
}

/// A result, in the one shape both revisions accept.
///
/// `resultType` is required by `2026-07-28` and unknown to `2025-11-25`, which ignores a field
/// it does not know — so writing it always is what lets one encoder serve both lanes.
fn result(id: &Value, result: Value) -> Value {
    let mut result = result;
    if let Some(object) = result.as_object_mut() {
        object.insert("resultType".to_string(), json!("complete"));
    }
    json!({ "jsonrpc": "2.0", "id": id, "result": result })
}

/// One tool, in the spelling MCP wants.
///
/// **The one rename is the only place the transport meets the schema**, so it is the only place
/// they can disagree — which is why the test compares the two lists field by field rather than
/// trusting this function.
fn tool_json(tool: &tools::Tool) -> Value {
    json!({
        "name": tool.name,
        "description": tool.description,
        "inputSchema": tool.input_schema,
    })
}

/// The name and version a client is told.
fn server_info() -> Value {
    json!({ "name": "azoth", "version": env!("CARGO_PKG_VERSION") })
}

/// Whether a message is a notification, which is one rule in one place because both transports
/// turn on it: it is answered with nothing on stdio and with `202 Accepted` on HTTP.
#[must_use]
pub fn is_notification(message: &Value) -> bool {
    message.get("id").is_none()
}

/// A JSON-RPC error, with the id echoed exactly as it arrived — a string id stays a string.
pub(crate) fn error_response(id: &Value, code: i64, message: &str, data: Option<Value>) -> Value {
    let mut error = json!({ "code": code, "message": message });
    if let Some(data) = data {
        error["data"] = data;
    }
    json!({ "jsonrpc": "2.0", "id": id, "error": error })
}

/// What a client is told about how to use the tools, in one paragraph.
const INSTRUCTIONS: &str = "\
Every tool call returns the whole flowsheet envelope: the document, the graph, the checker's \
diagnostics, every value's path, and the run. Reading a stream is reading the answer to the call \
that changed it, so there is no read tool. A call whose arguments are not a command answers \
`isError`; a call that landed and left the document refused answers `ok: false` with the \
diagnostics, which is a state to report rather than a failure; a run that failed is `run_error`.";
