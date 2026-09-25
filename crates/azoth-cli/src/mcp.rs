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
//! ignores a field it does not know. What differs is only which opening a client uses.
//!
//! **The transport is newline-delimited JSON with no embedded newlines.** `serde_json::to_string`
//! emits none, and every message is written and flushed as one line — a client blocks forever on a
//! partial one, which is the classic hand-rolled-framing failure and the reason the test reads a
//! line per request rather than writing everything and reading to the end.
//!
//! **Hand-rolled, and the cost is named.** Four things have to be right and only the last is cheap:
//! the framing above; error codes (`-32700` parse, `-32601` method, `-32602` tool, `-32022` the
//! protocol version, which is the specification's own number for it); the `input_schema` →
//! `inputSchema` rename; and the revision drift, which nothing here can detect. The class that
//! would make drift a *failing test* is a vendored `schema.json` for the revision plus a validator,
//! named in ROADMAP.md as not built rather than pretended away.

use std::io::{BufRead, Write};
use std::path::Path;

use azoth_core::{AzothError, Result};
use azoth_process::middleware::command::Command;
use azoth_process::middleware::session::Workspace;
use azoth_process::middleware::{envelope, tools};
use azoth_process::{UnitOpSpec, load_palette};
use serde_json::{Value, json};

/// The protocol revisions this server speaks, newest first.
///
/// **The newest is what a client that asks is told**, and both are what a version refusal lists.
/// Read from the specification rather than remembered: `2026-07-28` is current and handshake-free,
/// `2025-11-25` is the last handshake revision.
pub const PROTOCOL_VERSIONS: [&str; 2] = ["2026-07-28", "2025-11-25"];

/// JSON-RPC's own codes, plus the specification's one.
const PARSE_ERROR: i64 = -32700;
const METHOD_NOT_FOUND: i64 = -32601;
const UNKNOWN_TOOL: i64 = -32602;
const UNSUPPORTED_VERSION: i64 = -32022;

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
    let palette =
        load_palette(palette_dir).map_err(|error| AzothError::invalid_input("palette", error))?;
    let text = std::fs::read_to_string(flowsheet).map_err(|error| {
        AzothError::invalid_input("flowsheet", format!("{}: {error}", flowsheet.display()))
    })?;
    let workspace = Workspace::open(&text, palette.clone())
        .map_err(|error| AzothError::invalid_input("flowsheet", error.to_string()))?;

    let mut server = Server {
        workspace,
        palette,
        run: run_after_each_call,
    };

    for line in input.lines() {
        let line = line.map_err(|error| AzothError::invalid_input("stdin", error.to_string()))?;
        if line.trim().is_empty() {
            continue;
        }
        let response = match serde_json::from_str::<Value>(&line) {
            Ok(request) => server.handle(&request),
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

/// The session, and what this process was asked to do with it.
struct Server {
    workspace: Workspace,
    palette: Vec<UnitOpSpec>,
    run: bool,
}

impl Server {
    /// One request, as a response — or `None` for a notification.
    ///
    /// **A notification has no reply**, and answering one puts a message on the wire the client is
    /// not reading: the specification's own rule, and the reason this returns an `Option`.
    fn handle(&mut self, request: &Value) -> Option<Value> {
        let id = request.get("id").cloned();
        let method = request
            .get("method")
            .and_then(Value::as_str)
            .unwrap_or_default();

        // A request with no id is a notification, whatever its method: no reply, because a reply
        // would be a message the client is not reading.
        let id = id?;

        // The current revision carries the version on every request; a client that names one this
        // server does not speak is told which it could have used, per the specification.
        if let Some(version) = request
            .pointer("/params/_meta/io.modelcontextprotocol~1protocolVersion")
            .and_then(Value::as_str)
        {
            if !PROTOCOL_VERSIONS.contains(&version) {
                return Some(error_response(
                    &id,
                    UNSUPPORTED_VERSION,
                    &format!("`{version}` is not a protocol revision this server speaks"),
                    Some(json!({ "supported": PROTOCOL_VERSIONS })),
                ));
            }
        }

        match method {
            // The handshake a `2025-11-25` client opens with. Answered whatever opened before it,
            // because a second `initialize` is a client's business and not this server's.
            "initialize" => {
                let asked = request
                    .pointer("/params/protocolVersion")
                    .and_then(Value::as_str)
                    .unwrap_or(PROTOCOL_VERSIONS[0]);
                if !PROTOCOL_VERSIONS.contains(&asked) {
                    return Some(error_response(
                        &id,
                        UNSUPPORTED_VERSION,
                        &format!("`{asked}` is not a protocol revision this server speaks"),
                        Some(json!({ "supported": PROTOCOL_VERSIONS })),
                    ));
                }
                Some(self.result(
                    &id,
                    json!({
                        "protocolVersion": asked,
                        "capabilities": { "tools": { "listChanged": false } },
                        "serverInfo": server_info(),
                        "instructions": INSTRUCTIONS,
                    }),
                ))
            }
            // The current revision's mandatory discovery RPC.
            "server/discover" => Some(self.result(
                &id,
                json!({
                    "supportedVersions": PROTOCOL_VERSIONS,
                    "capabilities": { "tools": { "listChanged": false } },
                    "_meta": { "io.modelcontextprotocol/serverInfo": server_info() },
                    "instructions": INSTRUCTIONS,
                    "ttlMs": 3_600_000,
                    // A property of the software rather than of a document.
                    "cacheScope": "public",
                }),
            )),
            "ping" => Some(self.result(&id, json!({}))),
            "tools/list" => {
                let list: Vec<Value> = self.tools().iter().map(tool_json).collect();
                Some(self.result(
                    &id,
                    json!({
                        "tools": list,
                        "ttlMs": 3_600_000,
                        // A property of this build and this palette, so a client may cache it for
                        // the session but not across sessions.
                        "cacheScope": "private",
                    }),
                ))
            }
            "tools/call" => Some(self.call(&id, request)),
            other => Some(error_response(
                &id,
                METHOD_NOT_FOUND,
                &format!("`{other}` is not a method this server serves"),
                Some(json!({
                    "methods": ["initialize", "server/discover", "tools/list", "tools/call", "ping"]
                })),
            )),
        }
    }

    /// The tools, which are `middleware::tools`'s own list and nothing else.
    fn tools(&self) -> Vec<tools::Tool> {
        tools::tools(&self.palette)
    }

    /// One tool call: read the arguments as a command, apply it, answer with the envelope.
    ///
    /// **Three refusal states, kept apart because the middleware keeps them apart.** A call whose
    /// arguments are not a command is `isError`, because the model can correct it. A call that
    /// *landed* and left the document refused by the checker is **not** an error: a broken document
    /// is a state a caller shows, and `ok: false` is how it says so. A run that failed is not an
    /// error either — whether a substance exists is the databank's answer, and the envelope's
    /// `run_error` is where that arrives.
    fn call(&mut self, id: &Value, request: &Value) -> Value {
        let Some(name) = request.pointer("/params/name").and_then(Value::as_str) else {
            return error_response(id, UNKNOWN_TOOL, "a tool call must name a tool", None);
        };
        if !self.tools().iter().any(|tool| tool.name == name) {
            return error_response(id, UNKNOWN_TOOL, &format!("Unknown tool: {name}"), None);
        }

        // **The tool's name is the command's tag.** `tools.rs` writes it as the schema's `const`,
        // MCP carries the tool in `params.name`, and the command enum reads it from the object —
        // so this is the single point where the three spellings meet, and it is a copy of one
        // string rather than a fourth place it is written down.
        let mut arguments = request
            .pointer("/params/arguments")
            .cloned()
            .unwrap_or_else(|| json!({}));
        if let Some(object) = arguments.as_object_mut() {
            object.insert("command".to_string(), json!(name));
        }
        let command: Command = match serde_json::from_value(arguments) {
            Ok(command) => command,
            Err(error) => return self.refused(id, &error.to_string()),
        };
        if let Err(error) = self.workspace.apply(&command) {
            return self.refused(id, &error.to_string());
        }

        if self.run {
            // A run the edit made impossible is reported in the envelope rather than here, so the
            // error is deliberately dropped: it is the same fact, said in its own field.
            let _ = self.workspace.run();
        }

        match envelope::to_json(&self.workspace) {
            Ok(document) => {
                let structured: Value =
                    serde_json::from_str(&document).unwrap_or_else(|_| json!({}));
                self.result(
                    id,
                    json!({
                        "content": [ { "type": "text", "text": document } ],
                        // The same document twice: once for a model to read as text and once for a
                        // client to address field by field. One encoder wrote both.
                        "structuredContent": structured,
                        "isError": false,
                    }),
                )
            }
            Err(error) => self.refused(id, &error.to_string()),
        }
    }

    /// A call that could not take effect, as an `isError` result rather than a protocol error.
    ///
    /// The specification puts input-validation failures here rather than in a JSON-RPC error, so
    /// that a model can read the sentence and correct itself — and the sentence is the library's
    /// own, which names the offending field.
    fn refused(&self, id: &Value, reason: &str) -> Value {
        self.result(
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
    fn result(&self, id: &Value, result: Value) -> Value {
        let mut result = result;
        if let Some(object) = result.as_object_mut() {
            object.insert("resultType".to_string(), json!("complete"));
        }
        json!({ "jsonrpc": "2.0", "id": id, "result": result })
    }
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

/// A JSON-RPC error, with the id echoed exactly as it arrived — a string id stays a string.
fn error_response(id: &Value, code: i64, message: &str, data: Option<Value>) -> Value {
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
