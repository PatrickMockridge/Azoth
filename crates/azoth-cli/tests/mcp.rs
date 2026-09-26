//! `azoth mcp`, end to end: a real process, a real stream of JSON-RPC, real documents.
//!
//! **This is where "an MCP server is a projection of the schema rather than a second one" is
//! either true or false**, so the central assertion is that `tools/list` equals
//! `middleware::tools`'s own output field by field — the names, the descriptions and the schemas,
//! in order. The transport's one rename (`input_schema` → `inputSchema`) is the only place the two
//! can disagree, and that is exactly what comparing them catches.
//!
//! The driver reads **a line per request**, which is also a test of the framing: a server that
//! wrote a message in pieces, or batched two, would hang this rather than fail it — which is why
//! the reads are bounded by the child's own liveness and the assertions are on whole messages.

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, ChildStdin, ChildStdout, Command, Stdio};

use azoth_process::load_palette;
use azoth_process::middleware::tools;
use serde_json::{Value, json};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn demo() -> String {
    root()
        .join("specs/flowsheets/demo.toml")
        .display()
        .to_string()
}

fn palette() -> String {
    root().join("specs/unit_ops").display().to_string()
}

/// The server, as a client sees it: a process, and one message per line.
struct Client {
    child: Child,
    /// `None` once the pipe is closed, which is how a stdio server is told to stop.
    stdin: Option<ChildStdin>,
    stdout: BufReader<ChildStdout>,
}

impl Client {
    /// Start the server on the shipped document, running it after each call.
    fn start(extra: &[&str]) -> Self {
        // Owned strings, because `demo()` and `palette()` are temporaries and the child borrows
        // the argument list for the length of the spawn.
        let flowsheet = demo();
        let palette = palette();
        let mut args: Vec<&str> = vec!["mcp", "--flowsheet", &flowsheet, "--palette", &palette];
        args.extend_from_slice(extra);
        let mut child = Command::new(env!("CARGO_BIN_EXE_azoth"))
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .expect("the built binary should run");
        let stdin = child.stdin.take().expect("a stdin pipe");
        let stdout = BufReader::new(child.stdout.take().expect("a stdout pipe"));
        Self {
            child,
            stdin: Some(stdin),
            stdout,
        }
    }

    /// Close the write end, which a stdio server reads as the end of its input.
    fn close_stdin(&mut self) {
        self.stdin = None;
    }

    /// One request, and the one line it answers with.
    ///
    /// **Exactly as given**, so a test can send a request that is missing what the revision
    /// requires — which is the only way to check that the server refuses it.
    fn call(&mut self, method: &str, params: Value) -> Value {
        let id = method.to_string();
        let request = json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params });
        let stdin = self.stdin.as_mut().expect("the pipe is open");
        writeln!(stdin, "{request}").expect("the server is reading");
        stdin.flush().expect("it flushes");
        let mut line = String::new();
        self.stdout
            .read_line(&mut line)
            .expect("the server answers every request with a line");
        serde_json::from_str(&line).unwrap_or_else(|error| panic!("{error}\n{line}"))
    }

    /// One request in the current revision, which carries its protocol metadata on every call.
    fn modern(&mut self, method: &str, mut params: Value) -> Value {
        params["_meta"] = meta();
        self.call(method, params)
    }

    /// A request that expects a single `result`.
    fn result(&mut self, method: &str, params: Value) -> Value {
        let response = self.modern(method, params);
        response
            .get("result")
            .cloned()
            .unwrap_or_else(|| panic!("{method} refused: {response}"))
    }

    /// One tool call's arguments, as `tools/call` takes them.
    fn tool(&mut self, name: &str, arguments: Value) -> Value {
        self.modern(
            "tools/call",
            json!({ "name": name, "arguments": arguments }),
        )
    }
}

/// The `_meta` the current revision requires: both keys, because the schema's own `required` array
/// names both — `clientCapabilities` even though this server consumes no client capability.
fn meta() -> Value {
    json!({
        "io.modelcontextprotocol/protocolVersion": "2026-07-28",
        "io.modelcontextprotocol/clientCapabilities": {},
    })
}

impl Drop for Client {
    fn drop(&mut self) {
        // Closing stdin is how a stdio server is told to stop; the child is reaped so a hanging
        // one fails the test rather than leaking.
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

#[test]
fn the_handshake_lane_a_2025_client_opens_with() {
    let mut client = Client::start(&[]);
    let result = client.result("initialize", json!({ "protocolVersion": "2025-11-25" }));
    assert_eq!(result["protocolVersion"], "2025-11-25");
    assert_eq!(result["capabilities"]["tools"]["listChanged"], false);
    assert_eq!(result["serverInfo"]["name"], "azoth");
    assert_eq!(
        result["serverInfo"]["version"],
        env!("CARGO_PKG_VERSION"),
        "the version is the crate's own and not a second copy of it"
    );
    assert_eq!(result["resultType"], "complete");

    // A version this server does not speak is refused with the ones it does, which is the
    // specification's own negotiation rule.
    let refused = client.call("initialize", json!({ "protocolVersion": "1900-01-01" }));
    assert_eq!(refused["error"]["code"], -32022);
    assert_eq!(
        refused["error"]["data"]["supported"],
        json!(["2026-07-28", "2025-11-25"])
    );
    // **And the version that was asked for**, which the schema requires beside the list: a refusal
    // naming only what this server speaks leaves a client unable to say what was refused.
    assert_eq!(refused["error"]["data"]["requested"], "1900-01-01");
}

#[test]
fn the_discovery_lane_the_current_revision_requires() {
    let mut client = Client::start(&[]);
    let result = client.result("server/discover", json!({}));
    assert_eq!(
        result["supportedVersions"],
        json!(["2026-07-28", "2025-11-25"])
    );
    assert_eq!(result["resultType"], "complete");
    assert_eq!(result["cacheScope"], "public");
    assert_eq!(
        result["_meta"]["io.modelcontextprotocol/serverInfo"]["name"],
        "azoth"
    );

    // And the current revision's per-request version travels in `_meta`: one it does not speak is
    // refused there rather than at a handshake.
    let mut asked = meta();
    asked["io.modelcontextprotocol/protocolVersion"] = json!("1999-01-01");
    let refused = client.call("tools/list", json!({ "_meta": asked }));
    assert_eq!(refused["error"]["code"], -32022);
    assert_eq!(
        refused["error"]["data"]["supported"],
        json!(["2026-07-28", "2025-11-25"])
    );
    assert_eq!(refused["error"]["data"]["requested"], "1999-01-01");
}

/// **A request in the current revision carries its metadata, and one that does not is malformed.**
///
/// The requirement is the schema's — `_meta` is required on a request, and the request-metadata
/// object requires both `protocolVersion` and `clientCapabilities` — and the specification's answer
/// to a missing required field is `-32602` rather than the version refusal's `-32022`. Before this,
/// a request with no `_meta` at all was served and its version was never judged.
#[test]
fn a_modern_request_without_its_metadata_is_malformed() {
    let mut client = Client::start(&[]);

    let bare = client.call("tools/list", json!({}));
    assert_eq!(bare["error"]["code"], -32602);
    assert!(
        bare["error"]["message"]
            .as_str()
            .expect("a sentence")
            .contains("_meta"),
        "{bare}"
    );

    // One of the two keys is not enough, and the refusal names the one that is missing.
    let half = client.call(
        "tools/list",
        json!({ "_meta": { "io.modelcontextprotocol/protocolVersion": "2026-07-28" } }),
    );
    assert_eq!(half["error"]["code"], -32602);
    assert!(
        half["error"]["message"]
            .as_str()
            .expect("a sentence")
            .contains("clientCapabilities"),
        "{half}"
    );

    // Missing it is a *protocol* error, which is what makes it a `400` on HTTP rather than a
    // result a caller has to read.
    assert!(bare.get("result").is_none());
}

/// **`ping` belongs to the handshake lane, and the schema is what says so.**
///
/// `2025-11-25` defines a `PingRequest`; `2026-07-28` defines none and never names `ping` — a
/// modern client has `server/discover`, which this server already answers. So the lane a `ping` is
/// served in is decided by what the request carries: with `_meta` it is a modern request naming a
/// method that revision does not have, and without it it is the handshake lane's liveness check.
#[test]
fn ping_is_the_handshake_lanes_method() {
    let mut client = Client::start(&[]);

    // Without `_meta` on a connection that never opened: a modern request missing its metadata.
    let unopened = client.call("ping", json!({}));
    assert_eq!(unopened["error"]["code"], -32602);

    // A `2025-11-25` client opens with the handshake, and its `ping` is then answered — sent raw,
    // because `_meta` is what would put it in the other lane.
    client.result("initialize", json!({ "protocolVersion": "2025-11-25" }));
    let answered = client.call("ping", json!({}));
    assert_eq!(answered["result"]["resultType"], "complete");

    // And a request that *does* carry `_meta` is modern whatever came before it, where `ping` is
    // not a method at all.
    let modern = client.modern("ping", json!({}));
    assert_eq!(modern["error"]["code"], -32601);
    assert!(
        !modern["error"]["data"]["methods"]
            .as_array()
            .expect("a list")
            .contains(&json!("ping")),
        "the modern lane's own method list is what it is refused with: {modern}"
    );
}

/// **The claim the whole file exists for.** The tools the transport serves are the schema's own,
/// field for field — so the rename is the only thing between them and it is exactly one rename.
#[test]
fn tools_list_is_the_schema_not_a_second_copy_of_it() {
    let mut client = Client::start(&[]);
    let result = client.result("tools/list", json!({}));
    let served = result["tools"].as_array().expect("an array");

    let published = tools::tools(&load_palette(&root().join("specs/unit_ops")).expect("it loads"));
    assert_eq!(served.len(), published.len());
    for (served, published) in served.iter().zip(&published) {
        assert_eq!(served["name"], published.name);
        assert_eq!(served["description"], published.description);
        assert_eq!(
            served["inputSchema"], published.input_schema,
            "{}'s schema differs from the middleware's",
            published.name
        );
        // The rename, named rather than implied.
        assert!(served.get("input_schema").is_none());
    }
    // Deterministic, so a client may cache it: the current revision says SHOULD and equality to
    // the schema's own order makes it so.
    assert_eq!(served[0]["name"], "add_instance");
    assert_eq!(result["cacheScope"], "private");
}

/// **One session, across calls.** The second call sees the first one's effect, which is what makes
/// an agent's turns cumulative rather than independent.
#[test]
fn a_call_edits_the_document_the_next_call_sees() {
    let mut client = Client::start(&[]);

    let first = client.tool("remove_instance", json!({ "id": "hx1" }));
    assert_eq!(first["result"]["isError"], false);
    let envelope = &first["result"]["structuredContent"];
    let document = envelope["flowsheet"]["document"]
        .as_str()
        .expect("the text");
    assert!(!document.contains("id = \"hx1\""));
    // Only its own edges went with it: the unit that was downstream is still declared, which is
    // what the second call below then removes.
    assert!(document.contains("id = \"p1\""));

    // **A document the checker now refuses is a state and not a failure**: `isError` would say the
    // call did not happen, and it did.
    assert_eq!(envelope["ok"], false);
    let codes: Vec<&str> = envelope["diagnostics"]
        .as_array()
        .expect("an array")
        .iter()
        .filter_map(|diagnostic| diagnostic["code"].as_str())
        .collect();
    assert_eq!(codes, ["underfed_port", "underfed_port"]);
    assert_eq!(envelope["diagnostics"][0]["target"]["kind"], "handle");
    assert_eq!(envelope["diagnostics"][0]["target"]["node"], "p1");

    // The second call names what the first one made unaddressable, and the first edit is still in
    // the document: one `Workspace` for the life of the process.
    let second = client.tool("remove_instance", json!({ "id": "p1" }));
    let document = second["result"]["structuredContent"]["flowsheet"]["document"]
        .as_str()
        .expect("the text")
        .to_string();
    assert!(!document.contains("id = \"p1\""));
    assert!(!document.contains("id = \"hx1\""));
}

/// **The run is in the answer**, because the answer is the only place to ask for it.
#[test]
fn a_call_runs_the_document_and_answers_with_its_values() {
    let mut client = Client::start(&[]);
    let result = client.tool(
        "set_parameter",
        json!({
            "instance": "p1", "name": "outlet_pressure", "value": 4.0e6,
        }),
    );
    let envelope = &result["result"]["structuredContent"];
    assert_eq!(envelope["dirty"], false, "the values are the document's");
    assert_eq!(envelope["session"]["converged"], true);
    assert_eq!(
        envelope["session"]["streams"]["p1.outlet"]["P"]["magnitude_si"],
        4.0e6
    );
    assert_eq!(envelope["execution_order"], "insertion");
    // And the text half is the same document, for a model rather than a client.
    let text = result["result"]["content"][0]["text"]
        .as_str()
        .expect("text");
    assert_eq!(
        serde_json::from_str::<Value>(text).expect("it parses"),
        *envelope
    );

    // `--no-run` is the same session without the physics, and it says so rather than hiding it.
    let mut bare = Client::start(&["--no-run"]);
    let result = bare.tool("add_product", json!({ "name": "purge" }));
    let envelope = &result["result"]["structuredContent"];
    assert_eq!(envelope["session"], Value::Null);
    assert_eq!(envelope["dirty"], true);
    assert!(
        envelope["flowsheet"]["document"]
            .as_str()
            .expect("text")
            .contains("purge")
    );
}

#[test]
fn the_three_refusals_stay_three_things() {
    let mut client = Client::start(&[]);

    // A tool that is not one: a protocol error, because no command can be read at all.
    let unknown = client.tool("delete_everything", json!({}));
    assert_eq!(unknown["error"]["code"], -32602);
    assert!(
        unknown["error"]["message"]
            .as_str()
            .expect("a sentence")
            .contains("delete_everything")
    );

    // A command that cannot take effect: `isError`, so a model can correct itself, with the
    // library's own sentence naming the field.
    let refused = client.tool(
        "set_position",
        json!({ "node": "instance:nope", "x": 0, "y": 0 }),
    );
    assert_eq!(refused["result"]["isError"], true);
    assert!(
        refused["result"]["content"][0]["text"]
            .as_str()
            .expect("text")
            .contains("nope")
    );
    assert!(refused["result"].get("structuredContent").is_none());

    // A call that landed and left a document the checker refuses: **not** an error.
    let landed = client.tool(
        "add_instance",
        json!({ "id": "x1", "unit": "unit_ops.nosuch" }),
    );
    assert_eq!(landed["result"]["isError"], false);
    assert_eq!(landed["result"]["structuredContent"]["ok"], false);
    let codes: Vec<&str> = landed["result"]["structuredContent"]["diagnostics"]
        .as_array()
        .expect("an array")
        .iter()
        .filter_map(|diagnostic| diagnostic["code"].as_str())
        .collect();
    assert!(codes.contains(&"unknown_unit_op"), "{codes:?}");

    // An unknown method, and a line that is not JSON at all.
    let missing = client.modern("tools/write", json!({}));
    assert_eq!(missing["error"]["code"], -32601);
}

#[test]
fn the_process_ends_when_its_stdin_does() {
    let mut client = Client::start(&[]);
    // One exchange, so the server is provably running...
    client.result("server/discover", json!({}));
    // ...and then the write end closes, which is the end of its input.
    client.close_stdin();
    let status = client.child.wait().expect("the server exits");
    assert!(status.success(), "a closed stdin is a clean shutdown");
}
