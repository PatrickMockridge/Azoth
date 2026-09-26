//! `POST /mcp`, end to end: a real process, a real socket, real HTTP.
//!
//! The third of the three transport suites, and what it holds is what HTTP adds to the messages the
//! other two already test: the status a refusal is served with, the headers the specification
//! mirrors the body into, and the two things a modern-only endpoint owes a client speaking an older
//! revision. The messages themselves are `tests/mcp.rs`'s claim, and the calls are
//! `tests/serve.rs`'s — so nothing here re-asserts what those have shown.
//!
//! **One document, two doors**, is the exception, and it is here because it can only be here: an
//! edit through `/mcp` read back through `/rpc` is a claim about one `Session` behind two routes,
//! and neither of the other suites can make it.
//!
//! Every request is written and the response read to end of stream, which is what `Connection:
//! close` promises: a server that answered with a body and kept the connection open would hang this
//! rather than fail it.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpStream;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc;
use std::time::Duration;

use azoth_process::load_palette;
use azoth_process::middleware::tools;
use serde_json::{Value, json};

/// The revision every message in this file is written in, and the one the endpoint speaks.
const REV: &str = "2026-07-28";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The server, as a client sees it: a process, a port, and one request per connection.
struct Server {
    child: Child,
    port: u16,
}

impl Server {
    /// Start a server on a port the operating system picks, and read which one it got.
    fn start(origins: &[&str]) -> Self {
        let flowsheet = root()
            .join("specs/flowsheets/demo.toml")
            .display()
            .to_string();
        let palette = root().join("specs/unit_ops").display().to_string();
        let mut args = vec![
            "serve",
            "--flowsheet",
            &flowsheet,
            "--palette",
            &palette,
            "--port",
            "0",
        ];
        for origin in origins {
            args.push("--allow-origin");
            args.push(origin);
        }
        let mut child = Command::new(env!("CARGO_BIN_EXE_azoth"))
            .args(args)
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .expect("the built binary should run");

        // A reader thread and a deadline, because a server that failed to start would otherwise
        // block this test on a pipe nobody will write to.
        let stdout = child.stdout.take().expect("a stdout pipe");
        let (sender, receiver) = mpsc::channel();
        std::thread::spawn(move || {
            for line in BufReader::new(stdout).lines().map_while(Result::ok) {
                if sender.send(line).is_err() {
                    return;
                }
            }
        });
        let mut port = None;
        while let Ok(line) = receiver.recv_timeout(Duration::from_secs(30)) {
            if let Some(rest) = line.strip_prefix("azoth: serving http://") {
                let address = rest.split('/').next().unwrap_or_default();
                port = address.rsplit(':').next().and_then(|p| p.parse().ok());
                break;
            }
        }
        Self {
            child,
            port: port.expect("the server prints the address it bound"),
        }
    }

    /// One request, and the status, the response head and the body it answers with.
    fn request(
        &self,
        method: &str,
        path: &str,
        headers: &[(String, String)],
        body: &str,
    ) -> (u16, String, Value) {
        let mut stream = TcpStream::connect(("127.0.0.1", self.port)).expect("it is listening");
        let mut head = format!("{method} {path} HTTP/1.1\r\nHost: 127.0.0.1\r\n");
        for (key, value) in headers {
            head.push_str(&format!("{key}: {value}\r\n"));
        }
        head.push_str(&format!("Content-Length: {}\r\n\r\n", body.len()));
        stream
            .write_all(head.as_bytes())
            .and_then(|()| stream.write_all(body.as_bytes()))
            .expect("it writes");
        stream.flush().expect("it flushes");

        let mut response = String::new();
        stream
            .read_to_string(&mut response)
            .expect("the server closes the connection after one answer");
        let (head, body) = response
            .split_once("\r\n\r\n")
            .expect("a response has a head and a body");
        let status: u16 = head
            .split_whitespace()
            .nth(1)
            .and_then(|code| code.parse().ok())
            .expect("a status line");
        let parsed = if body.trim().is_empty() {
            Value::Null
        } else {
            serde_json::from_str(body).unwrap_or_else(|error| panic!("{error}\n{body}"))
        };
        (status, head.to_string(), parsed)
    }

    /// One message, with the headers a conforming client derives from it.
    fn send(&self, message: &Value) -> (u16, Value) {
        let (status, _, body) =
            self.request("POST", "/mcp", &headers_of(message), &message.to_string());
        (status, body)
    }

    /// One message with the headers replaced, for the cases where a client's headers lie.
    fn send_with(&self, message: &Value, headers: &[(&str, &str)]) -> (u16, Value) {
        let headers: Vec<(String, String)> = headers
            .iter()
            .map(|(key, value)| ((*key).to_string(), (*value).to_string()))
            .collect();
        let (status, _, body) = self.request("POST", "/mcp", &headers, &message.to_string());
        (status, body)
    }

    /// One request the endpoint should refuse before it reads a message at all.
    fn send_raw(&self, method: &str, headers: &[(&str, &str)]) -> (u16, String, Value) {
        let headers: Vec<(String, String)> = headers
            .iter()
            .map(|(key, value)| ((*key).to_string(), (*value).to_string()))
            .collect();
        self.request(method, "/mcp", &headers, "")
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

/// One message in the current revision, which carries its metadata in `_meta`.
fn message(id: i64, method: &str, mut params: Value) -> Value {
    params["_meta"] = json!({
        "io.modelcontextprotocol/protocolVersion": REV,
        "io.modelcontextprotocol/clientCapabilities": {},
    });
    json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params })
}

/// The headers a conforming client sends for a message, derived from the body exactly as the
/// specification says a client derives them — which is what makes `send_with` a *lie* rather than a
/// second way to send one.
fn headers_of(message: &Value) -> Vec<(String, String)> {
    let mut headers = vec![
        ("MCP-Protocol-Version".to_string(), REV.to_string()),
        (
            "Mcp-Method".to_string(),
            message["method"].as_str().unwrap_or_default().to_string(),
        ),
    ];
    if let Some(name) = message.pointer("/params/name").and_then(Value::as_str) {
        headers.push(("Mcp-Name".to_string(), name.to_string()));
    }
    headers
}

/// The error code in a JSON-RPC error body, which every refusal below names explicitly.
fn code(body: &Value) -> i64 {
    body.pointer("/error/code")
        .and_then(Value::as_i64)
        .unwrap_or_else(|| panic!("no error code in {body}"))
}

/// **The same three methods the stdio door serves, over HTTP.**
#[test]
fn the_catalogue_a_discovery_and_a_call() {
    let server = Server::start(&[]);

    let (status, discover) = server.send(&message(1, "server/discover", json!({})));
    assert_eq!(status, 200);
    assert_eq!(discover["result"]["resultType"], "complete");
    assert_eq!(discover["result"]["cacheScope"], "public");

    let (status, list) = server.send(&message(2, "tools/list", json!({})));
    assert_eq!(status, 200);
    let served = list["result"]["tools"].as_array().expect("an array");
    let published = tools::tools(&load_palette(&root().join("specs/unit_ops")).expect("it loads"));
    assert_eq!(served.len(), published.len());
    assert_eq!(served[0]["name"], "add_instance");

    let (status, called) = server.send(&message(
        3,
        "tools/call",
        json!({ "name": "set_parameter", "arguments": {
            "instance": "p1", "name": "outlet_pressure", "value": 4.0e6,
        }}),
    ));
    assert_eq!(status, 200);
    let envelope = &called["result"]["structuredContent"];
    assert_eq!(envelope["dirty"], false);
    assert_eq!(
        envelope["session"]["streams"]["p1.outlet"]["P"]["magnitude_si"],
        4.0e6
    );
}

/// **One document behind two doors**, which is why the endpoint is on `azoth serve` at all.
#[test]
fn an_edit_through_mcp_is_read_back_through_rpc() {
    let server = Server::start(&[]);

    let (status, _) = server.send(&message(
        1,
        "tools/call",
        json!({ "name": "remove_instance", "arguments": { "id": "hx1" } }),
    ));
    assert_eq!(status, 200);

    // The editor's door, and the same `Session`: the unit the agent removed is gone from the
    // document this reads.
    let (status, head, envelope) = server.request("POST", "/rpc", &[], "{}");
    assert_eq!(status, 200, "{envelope}");
    assert!(head.starts_with("HTTP/1.1 200"), "{head}");
    let document = envelope["flowsheet"]["document"]
        .as_str()
        .expect("the text");
    assert!(!document.contains("id = \"hx1\""));
    assert!(document.contains("id = \"p1\""));
}

/// **The header check, which is the one rule this transport exists to get right.**
///
/// An intermediary routes on `Mcp-Method` and this process executes the body's, so the two must be
/// the same message — and a disagreement is `-32020` with a `400`, which the schema states outright:
/// *"For HTTP, the response status code MUST be `400 Bad Request`."*
#[test]
fn a_header_that_disagrees_with_the_body_is_refused() {
    let server = Server::start(&[]);
    let list = message(1, "tools/list", json!({}));

    // The version header, missing and disagreeing.
    let (status, body) = server.send_with(
        &list,
        &[
            ("Mcp-Method", "tools/list"),
            ("MCP-Protocol-Version", "1999-01-01"),
        ],
    );
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32020);

    let (status, body) = server.send_with(&list, &[("Mcp-Method", "tools/list")]);
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32020);
    assert!(
        body["error"]["message"]
            .as_str()
            .expect("a sentence")
            .contains("MCP-Protocol-Version"),
        "{body}"
    );

    // The method header, missing and disagreeing.
    let (status, body) = server.send_with(&list, &[("MCP-Protocol-Version", REV)]);
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32020);

    let (status, body) = server.send_with(
        &list,
        &[("MCP-Protocol-Version", REV), ("Mcp-Method", "tools/call")],
    );
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32020);
    assert!(
        body["error"]["message"]
            .as_str()
            .expect("a sentence")
            .contains("tools/call"),
        "{body}"
    );

    // And the name, for the one request whose name the specification requires a header for.
    let call = message(
        2,
        "tools/call",
        json!({ "name": "remove_instance", "arguments": { "id": "hx1" } }),
    );
    let (status, body) = server.send_with(
        &call,
        &[
            ("MCP-Protocol-Version", REV),
            ("Mcp-Method", "tools/call"),
            ("Mcp-Name", "add_instance"),
        ],
    );
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32020);
    assert!(
        body["error"]["message"]
            .as_str()
            .expect("a sentence")
            .contains("add_instance"),
        "{body}"
    );
}

/// **An encoded header is decoded before it is compared**, which is the specification's rule for a
/// server that inspects these values — a tool name outside the header-safe set travels Base64, and
/// comparing the encoded text against the body would refuse a conforming client.
#[test]
fn an_encoded_name_header_is_decoded_before_it_is_compared() {
    let server = Server::start(&[]);
    let call = message(
        1,
        "tools/call",
        json!({ "name": "remove_instance", "arguments": { "id": "hx1" } }),
    );

    // `remove_instance`, encoded: the name is ASCII, so a real client would send it plain — but the
    // specification's rule is that a server decodes what it is given, so this must be *served*.
    let (status, body) = server.send_with(
        &call,
        &[
            ("MCP-Protocol-Version", REV),
            ("Mcp-Method", "tools/call"),
            ("Mcp-Name", "=?base64?cmVtb3ZlX2luc3RhbmNl?="),
        ],
    );
    assert_eq!(status, 200, "{body}");
    assert_eq!(body["result"]["isError"], false);

    // And an encoded name that decodes to something else is the same disagreement as any other.
    let (status, body) = server.send_with(
        &call,
        &[
            ("MCP-Protocol-Version", REV),
            ("Mcp-Method", "tools/call"),
            ("Mcp-Name", "=?base64?YWRkX2luc3RhbmNl?="),
        ],
    );
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32020);
}

/// **The metadata the revision requires, and the status the specification pins for its absence.**
#[test]
fn a_request_without_its_metadata_is_a_bad_request() {
    let server = Server::start(&[]);

    // No `_meta` at all, and a version header: the header has nothing in the body to agree with, so
    // this is the *header* refusal rather than the malformed-request one. Both are `400` — which is
    // what the specification pins — and the order the transport checks them in is what decides
    // which sentence a client gets.
    let bare = json!({ "jsonrpc": "2.0", "id": 1, "method": "tools/list", "params": {} });
    let (status, body) = server.send(&bare);
    assert_eq!(
        status, 400,
        "the specification says HTTP 400 for this: {body}"
    );
    assert_eq!(code(&body), -32020);
    assert!(
        body["error"]["message"]
            .as_str()
            .expect("a sentence")
            .contains("names nothing"),
        "{body}"
    );

    // Half of it: the schema requires both keys together, and here the version header *can* agree
    // with the body — so the disagreement is the body's own and it is the malformed-request code.
    let half = json!({ "jsonrpc": "2.0", "id": 2, "method": "tools/list",
        "params": { "_meta": { "io.modelcontextprotocol/protocolVersion": REV } } });
    let (status, body) = server.send(&half);
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32602);
    assert!(
        body["error"]["message"]
            .as_str()
            .expect("a sentence")
            .contains("clientCapabilities"),
        "{body}"
    );

    // An unsupported version is the other refusal, and it names both the versions this server
    // speaks and the one that was asked for.
    let mut wrong = message(3, "tools/list", json!({}));
    wrong["params"]["_meta"]["io.modelcontextprotocol/protocolVersion"] = json!("1999-01-01");
    let (status, body) = server.send_with(
        &wrong,
        &[
            ("MCP-Protocol-Version", "1999-01-01"),
            ("Mcp-Method", "tools/list"),
        ],
    );
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32022);
    assert_eq!(body["error"]["data"]["requested"], "1999-01-01");
    assert_eq!(
        body["error"]["data"]["supported"],
        json!([REV, "2025-11-25"])
    );
}

/// **The handshake this endpoint does not speak, refused with the revision it does.**
///
/// The specification asks a modern-only server to name its versions in any error it returns to an
/// `initialize`, because a legacy client has no fall-forward. **The list is this transport's**, so
/// it names one version and not the two the process speaks — a client told `2025-11-25` here would
/// retry with a revision this door refuses again.
#[test]
fn the_handshake_is_refused_with_where_that_lane_lives() {
    let server = Server::start(&[]);
    let initialize = json!({ "jsonrpc": "2.0", "id": 1, "method": "initialize",
        "params": { "protocolVersion": "2025-11-25" } });

    // **No headers at all**, because that is what a legacy client sends: its revision never defined
    // `MCP-Protocol-Version` or `Mcp-Method`, and its version travelled in the handshake's own
    // `params`. A server that judged this request by the header rule would refuse a conforming
    // client for a header it had no way to send.
    let (status, head, body) = server.request("POST", "/mcp", &[], &initialize.to_string());
    assert_eq!(status, 400, "{body}");
    assert_eq!(code(&body), -32022);
    assert_eq!(body["error"]["data"]["supported"], json!([REV]));
    assert_eq!(body["error"]["data"]["requested"], "2025-11-25");
    let message = body["error"]["message"].as_str().expect("a sentence");
    assert!(message.contains("azoth mcp"), "{message}");
    assert!(message.contains("stdio"), "{message}");
    assert!(head.starts_with("HTTP/1.1 400"), "{head}");
}

/// **An unknown method is `404`, and every other refusal is `400`** — the one status the
/// specification maps differently from the rest, and the transport's own mapping.
#[test]
fn a_method_this_revision_does_not_have_is_not_found() {
    let server = Server::start(&[]);

    let (status, body) = server.send(&message(1, "tools/write", json!({})));
    assert_eq!(status, 404, "{body}");
    assert_eq!(code(&body), -32601);

    // `ping` is the handshake lane's, which the schema says by not defining it.
    let (status, body) = server.send(&message(2, "ping", json!({})));
    assert_eq!(status, 404);
    assert_eq!(code(&body), -32601);

    // A tool that is not one is *not* a 404: the method is served, and `-32602` is what the
    // protocol says for an unknown tool.
    let (status, body) = server.send(&message(3, "tools/call", json!({ "name": "nope" })));
    assert_eq!(status, 400);
    assert_eq!(code(&body), -32602);
}

/// **The two methods the revision removed**, answered `405` because that is what the specification
/// tells a server in this position to answer them.
#[test]
fn the_stream_and_the_session_teardown_are_gone() {
    let server = Server::start(&[]);

    for method in ["GET", "DELETE"] {
        let (status, head, _) = server.send_raw(method, &[("MCP-Protocol-Version", REV)]);
        assert_eq!(status, 405, "{method}: {head}");
        assert!(
            head.starts_with("HTTP/1.1 405 Method Not Allowed"),
            "{head}"
        );
    }

    // And the headers of the revisions that had them are ignored rather than honoured: this server
    // mints no session and resumes no stream, so a client sending either is served normally.
    let (status, body) = server.send_with(
        &message(1, "tools/list", json!({})),
        &[
            ("MCP-Protocol-Version", REV),
            ("Mcp-Method", "tools/list"),
            ("Mcp-Session-Id", "deadbeef"),
            ("Last-Event-ID", "7"),
        ],
    );
    assert_eq!(status, 200, "{body}");
    assert!(body["result"]["tools"].is_array());
}

/// **A notification is accepted and not answered**, which on HTTP is a status rather than a
/// silence.
#[test]
fn a_notification_is_accepted_with_nothing_in_it() {
    let server = Server::start(&[]);
    let notification = json!({ "jsonrpc": "2.0", "method": "notifications/progress",
        "params": { "progressToken": "t", "progress": 1 } });

    // **No headers and no `_meta`**: the specification defines neither for a notification POST, so
    // checking them here would refuse a message it says to accept.
    let (status, head, body) = server.request("POST", "/mcp", &[], &notification.to_string());
    assert_eq!(status, 202, "{head}");
    assert!(head.starts_with("HTTP/1.1 202 Accepted"), "{head}");
    assert_eq!(body, Value::Null, "and nothing is written to read");
    assert!(
        !head.contains("Content-Type"),
        "a response with no body has no content to describe: {head}"
    );
}

/// **A JSON-only endpoint refuses honestly**, rather than answering with a body the client has said
/// it cannot read.
#[test]
fn a_client_that_cannot_read_json_is_told_so() {
    let server = Server::start(&[]);
    let list = message(1, "tools/list", json!({}));
    let headers = headers_of(&list);

    let mut refused = headers.clone();
    refused.push(("Accept".to_string(), "text/event-stream".to_string()));
    let (status, _, body) = server.request("POST", "/mcp", &refused, &list.to_string());
    assert_eq!(status, 406, "{body}");

    // Both, which is what the specification requires a client to send, and a wildcard, are served.
    for accept in [
        "application/json, text/event-stream",
        "*/*",
        "application/json",
    ] {
        let mut served = headers.clone();
        served.push(("Accept".to_string(), accept.to_string()));
        let (status, _, body) = server.request("POST", "/mcp", &served, &list.to_string());
        assert_eq!(status, 200, "`{accept}`: {body}");
    }
}

/// **The origin policy is the server's, and this door makes it a refusal.**
///
/// `/rpc` withholds the CORS header from a disallowed origin and processes the call anyway, which is
/// what a browser reads as "no". The specification requires *this* transport to answer `403`, and
/// the two behaviours are deliberately different.
#[test]
fn a_disallowed_origin_is_forbidden_here() {
    let allowed = "http://localhost:5173";
    let server = Server::start(&[allowed]);
    let list = message(1, "tools/list", json!({}));
    let mut headers = headers_of(&list);
    headers.push(("Origin".to_string(), "http://evil.example".to_string()));
    let (status, head, _) = server.request("POST", "/mcp", &headers, &list.to_string());
    assert_eq!(status, 403, "{head}");
    assert!(
        !head.contains("Access-Control-Allow-Origin"),
        "a refusal a browser must read as one: {head}"
    );

    let mut headers = headers_of(&list);
    headers.push(("Origin".to_string(), allowed.to_string()));
    let (status, head, body) = server.request("POST", "/mcp", &headers, &list.to_string());
    assert_eq!(status, 200, "{body}");
    assert!(
        head.contains(&format!("Access-Control-Allow-Origin: {allowed}")),
        "{head}"
    );

    // A client with no `Origin` at all is served, which is what makes `curl` work and is the half
    // of this policy the editor's door shares.
    let (status, _) = server.send(&list);
    assert_eq!(status, 200);
}

/// **A body that is not a message is a parse error**, with the id the specification gives it: none.
#[test]
fn a_body_that_is_not_json_has_no_id_to_answer_to() {
    let server = Server::start(&[]);

    let (status, _, body) = server.request("POST", "/mcp", &[], "{not json");
    assert_eq!(status, 400, "{body}");
    assert_eq!(code(&body), -32700);
    assert_eq!(body["id"], Value::Null);
}
