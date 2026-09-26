//! `azoth serve`, end to end: a real process, a real socket, real HTTP.
//!
//! The sibling of `tests/mcp.rs`, and what it holds is the same claim on the other transport:
//! **the calls are the session's and the answer is the envelope**, so a client that speaks one
//! transport speaks the other. What is *only* this file's is the parts HTTP owns - the route, the
//! status codes, the origin policy and the framing - and those are what the assertions below are
//! about, because a second transport is where a second set of rules would hide.
//!
//! Every request is written and the response read to end of stream, which is what
//! `Connection: close` promises: a server that answered with a body and kept the connection open
//! would hang this rather than fail it.

use std::io::{BufRead, BufReader, Read, Write};
use std::net::TcpStream;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::mpsc;
use std::time::Duration;

use serde_json::{Value, json};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The server, as a client sees it: a process, a port, and one request per connection.
struct Server {
    child: Child,
    port: u16,
    origins: Vec<String>,
}

impl Server {
    /// Start a server on a port the operating system picks, and read which one it got.
    ///
    /// **The printed line is the contract**, which is why the test reads it rather than picking a
    /// port: a fixed port makes two test binaries collide, and `--port 0` means only the server
    /// knows what it bound.
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
            origins: origins.iter().map(|o| (*o).to_string()).collect(),
        }
    }

    /// One request, and the status and body it answers with.
    fn request(&self, method: &str, body: Option<&str>, origin: Option<&str>) -> (u16, Value) {
        let mut stream = TcpStream::connect(("127.0.0.1", self.port)).expect("it is listening");
        let body = body.unwrap_or_default();
        let origin = origin
            .map(|origin| format!("Origin: {origin}\r\n"))
            .unwrap_or_default();
        let request = format!(
            "{method} /rpc HTTP/1.1\r\nHost: 127.0.0.1\r\n{origin}\
             Content-Type: application/json\r\nContent-Length: {}\r\n\r\n{body}",
            body.len()
        );
        stream.write_all(request.as_bytes()).expect("it writes");
        stream.flush().expect("it flushes");

        let mut response = String::new();
        stream
            .read_to_string(&mut response)
            .expect("the server closes the connection after one answer");
        let (head, body) = response
            .split_once("\r\n\r\n")
            .expect("a response has a body");
        let status: u16 = head
            .split_whitespace()
            .nth(1)
            .and_then(|code| code.parse().ok())
            .expect("a status line");
        let parsed = if body.trim().is_empty() {
            Value::Null
        } else {
            serde_json::from_str(body).expect("the body is JSON")
        };
        (status, parsed)
    }

    /// The headers of one response, for the origin policy.
    fn headers(&self, method: &str, origin: Option<&str>) -> String {
        let mut stream = TcpStream::connect(("127.0.0.1", self.port)).expect("it is listening");
        let origin = origin
            .map(|origin| format!("Origin: {origin}\r\n"))
            .unwrap_or_default();
        let request = format!(
            "{method} /rpc HTTP/1.1\r\nHost: 127.0.0.1\r\n{origin}\
             Access-Control-Request-Method: POST\r\nContent-Length: 0\r\n\r\n"
        );
        stream.write_all(request.as_bytes()).expect("it writes");
        let mut response = String::new();
        stream.read_to_string(&mut response).expect("it answers");
        response
            .split_once("\r\n\r\n")
            .map(|(head, _)| head.to_string())
            .unwrap_or_default()
    }

    /// One call, which every test below makes by name.
    fn call(&self, body: Value) -> (u16, Value) {
        self.request("POST", Some(&body.to_string()), None)
    }

    /// A call that is expected to succeed, returning the envelope.
    fn envelope(&self, body: Value) -> Value {
        let (status, answer) = self.call(body);
        assert_eq!(status, 200, "the call is served: {answer}");
        answer
    }
}

impl Drop for Server {
    fn drop(&mut self) {
        let _ = self.child.kill();
        let _ = self.child.wait();
    }
}

#[test]
fn the_envelope_is_the_answer_to_a_look() {
    let server = Server::start(&[]);
    assert!(server.origins.is_empty());
    let envelope = server.envelope(json!({}));
    assert_eq!(envelope["ok"], true);
    assert_eq!(envelope["flowsheet"]["id"], "flowsheets.demo");
    assert_eq!(
        envelope["flowsheet"]["graph"]["nodes"]
            .as_array()
            .map(Vec::len),
        Some(6)
    );
    // A look does not run: the document has not been solved, and saying so is what `dirty` is.
    assert_eq!(envelope["dirty"], true);
    assert_eq!(envelope["session"], Value::Null);
}

/// **One session, across requests.** The second call sees the first one's effect, which is what
/// makes a hosted document a document rather than a sequence of unrelated edits.
#[test]
fn a_call_edits_the_document_the_next_call_sees() {
    let server = Server::start(&[]);

    let first = server.envelope(json!({
        "command": {"command": "remove_instance", "id": "hx1"},
    }));
    let document = first["flowsheet"]["document"].as_str().expect("the text");
    assert!(!document.contains("id = \"hx1\""));
    // **A document the checker now refuses is a state and not a failure**: `ok: false` with the
    // diagnostics, because the call happened.
    assert_eq!(first["ok"], false);
    let codes: Vec<&str> = first["diagnostics"]
        .as_array()
        .expect("an array")
        .iter()
        .filter_map(|diagnostic| diagnostic["code"].as_str())
        .collect();
    assert_eq!(codes, ["underfed_port", "underfed_port"]);

    let second = server.envelope(json!({
        "command": {"command": "remove_instance", "id": "p1"},
    }));
    let document = second["flowsheet"]["document"].as_str().expect("the text");
    assert!(!document.contains("id = \"p1\""));
    assert!(
        !document.contains("id = \"hx1\""),
        "the first edit is still there"
    );
}

/// **The palette is the one thing a hosted editor has nowhere else to ask for.** A form per unit op
/// is on no envelope, so a client that draws a widget cannot read it off the answer to a call —
/// which is why this is a request of its own, and why it is a *read*: nothing here is handed a
/// document.
#[test]
fn the_catalogue_is_the_palette_the_process_loaded() {
    let server = Server::start(&[]);

    let (status, catalogue) = server.call(json!({ "catalogue": false }));
    assert_eq!(status, 200);
    assert_eq!(
        catalogue["unit_ops"].as_array().map(Vec::len),
        Some(29),
        "the shipped palette, entry by entry"
    );
    // A form carries the parameters a unit-op window is drawn from, which is the whole point of
    // asking: without them a selected node has no fields.
    let pump = catalogue["unit_ops"]
        .as_array()
        .expect("an array")
        .iter()
        .find(|form| form["id"] == "unit_ops.pump")
        .expect("the palette has a pump");
    assert!(
        !pump["parameters"]
            .as_array()
            .expect("parameters")
            .is_empty(),
        "{pump}"
    );
    // The tools are the agent's schema, and they are asked for rather than always carried.
    assert_eq!(catalogue["tools"], Value::Null);

    let (status, catalogue) = server.call(json!({ "catalogue": true }));
    assert_eq!(status, 200);
    assert_eq!(
        catalogue["tools"].as_array().map(Vec::len),
        Some(15),
        "the fifteen commands, which `/mcp` serves the same way"
    );

    // **A body that asks for two things has two answers**, so it is a request this server cannot
    // read rather than one it serves by preferring a field.
    let (status, answer) = server.call(json!({ "catalogue": false, "command": {} }));
    assert_eq!(status, 400);
    assert!(
        answer["error"]
            .as_str()
            .expect("a sentence")
            .contains("nothing else"),
        "{answer}"
    );

    let (status, answer) = server.call(json!({ "catalogue": "yes" }));
    assert_eq!(status, 400);
    assert!(
        answer["error"]
            .as_str()
            .expect("a sentence")
            .contains("true or false"),
        "{answer}"
    );
}

#[test]
fn an_edit_runs_the_document_unless_the_request_says_not_to() {
    let server = Server::start(&[]);
    let edit = json!({
        "command": {"command": "set_parameter", "instance": "p1",
                    "name": "outlet_pressure", "value": 4.0e6},
    });

    let mut without = edit.clone();
    without["run"] = json!(false);
    let envelope = server.envelope(without);
    assert_eq!(
        envelope["session"],
        Value::Null,
        "the edit did not pay the physics"
    );
    assert_eq!(envelope["dirty"], true);

    let envelope = server.envelope(edit);
    assert_eq!(envelope["dirty"], false);
    assert_eq!(
        envelope["session"]["streams"]["p1.outlet"]["P"]["magnitude_si"],
        4.0e6
    );
}

#[test]
fn a_look_runs_only_when_it_asks_to() {
    let server = Server::start(&[]);
    // Nothing has run yet, and a look says so rather than running silently.
    assert_eq!(server.envelope(json!({}))["session"], Value::Null);
    // Asking for one is a different request with the same shape.
    let envelope = server.envelope(json!({ "run": true }));
    assert_eq!(envelope["session"]["converged"], true);
    assert_eq!(envelope["session"]["iterations"], 2);
    assert_eq!(envelope["dirty"], false);
}

#[test]
fn the_order_is_settable_and_is_not_a_change_to_the_document() {
    let server = Server::start(&[]);
    let envelope = server.envelope(json!({ "order": "topological" }));
    assert_eq!(envelope["execution_order"], "topological");
    assert_eq!(envelope["dirty"], true, "an order is not an edit");
    assert_eq!(envelope["ok"], true, "and it is not a refusal either");

    // A name that is neither is a request this server cannot serve.
    let (status, answer) = server.call(json!({ "order": "kahn" }));
    assert_eq!(status, 400);
    assert!(
        answer["error"]
            .as_str()
            .expect("a sentence")
            .contains("kahn"),
        "{answer}"
    );
}

/// The same three states the MCP transport keeps apart, in HTTP's own terms.
#[test]
fn the_three_refusals_stay_three_things() {
    let server = Server::start(&[]);

    // A tool that is not one, a command that cannot take effect: the *call* is wrong, so the
    // client is told to fix it and the body carries the library's own sentence.
    let (status, answer) = server.call(json!({ "command": {"command": "delete_everything"} }));
    assert_eq!(status, 400);
    assert!(
        answer["error"]
            .as_str()
            .expect("a sentence")
            .contains("Unknown tool")
    );

    let (status, answer) = server.call(json!({
        "command": {"command": "set_position", "node": "instance:nope", "x": 0, "y": 0},
    }));
    assert_eq!(status, 400);
    assert!(
        answer["error"]
            .as_str()
            .expect("a sentence")
            .contains("nope")
    );

    // A call that *landed* and left a document the checker refuses: served, with `ok: false`.
    let (status, envelope) = server.call(json!({
        "command": {"command": "add_instance", "id": "x1", "unit": "unit_ops.nosuch"},
    }));
    assert_eq!(status, 200);
    assert_eq!(envelope["ok"], false);
    let codes: Vec<&str> = envelope["diagnostics"]
        .as_array()
        .expect("an array")
        .iter()
        .filter_map(|diagnostic| diagnostic["code"].as_str())
        .collect();
    assert!(codes.contains(&"unknown_unit_op"), "{codes:?}");
}

#[test]
fn a_request_this_server_cannot_read_is_a_bad_request_rather_than_a_call() {
    let server = Server::start(&[]);

    // **A body that parses and is the wrong shape**, which is a malformed *request* and not a call
    // that failed: nothing here names a document, so there is nothing to answer with but the
    // reason. The distinction the three refusal states turn on is whether a call *happened*.
    let (status, answer) = server.call(json!("not an object"));
    assert_eq!(status, 400);
    assert!(
        answer["error"]
            .as_str()
            .expect("a sentence")
            .contains("not an object"),
        "{answer}"
    );

    // And a body that is not JSON at all, written raw because `json!` cannot produce one.
    let (status, answer) = server.request("POST", Some("{not json"), None);
    assert_eq!(status, 400);
    assert!(
        answer["error"]
            .as_str()
            .expect("a sentence")
            .contains("not JSON"),
        "{answer}"
    );

    // An edit with no `command` field names no tool.
    let (status, answer) = server.call(json!({ "command": {"id": "hx1"} }));
    assert_eq!(status, 400);
    assert!(answer["error"].as_str().is_some(), "{answer}");

    // A route that is not the one route.
    let (status, _) = server.request("GET", None, None);
    assert_eq!(status, 404);
}

/// **A response without the CORS header is one a page cannot read**, so an origin a server was not
/// told to allow is served by nothing in a browser - and `curl` is unaffected, which is the point
/// of the posture rather than an accident of it.
#[test]
fn only_an_allowed_origin_is_given_the_cors_header() {
    let allowed = "http://localhost:5173";
    let server = Server::start(&[allowed]);
    assert_eq!(server.origins, [allowed]);

    let head = server.headers("OPTIONS", Some(allowed));
    assert!(head.starts_with("HTTP/1.1 204"), "{head}");
    assert!(
        head.contains(&format!("Access-Control-Allow-Origin: {allowed}")),
        "{head}"
    );

    // A page on another origin gets the refusal *without* the header, which is what the browser
    // reads as "no".
    let head = server.headers("OPTIONS", Some("http://evil.example"));
    assert!(head.starts_with("HTTP/1.1 403"), "{head}");
    assert!(!head.contains("Access-Control-Allow-Origin"), "{head}");

    // And a POST carries it too, so the answer to a real call is readable by the allowed page.
    let head = server.headers("POST", Some(allowed));
    assert!(
        head.contains(&format!("Access-Control-Allow-Origin: {allowed}")),
        "{head}"
    );
    assert!(
        head.contains("Connection: close"),
        "every response says so: {head}"
    );
}
