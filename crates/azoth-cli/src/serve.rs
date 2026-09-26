//! `azoth serve` — the same session, over HTTP.
//!
//! **A transport and not a second surface.** The calls are [`crate::session::Session`]'s, the answer
//! is the envelope, and this file is a request parser and a response writer: it decides nothing
//! about a flowsheet. It serves two doors — the editor's `/rpc` and the agent's `/mcp` — over one
//! session, and `crate::mcp_http` owns the second door's protocol.
//!
//! **Why a hosted session exists at all**, given that the library is the backend and the browser
//! runs the kernels itself as wasm: a *hosted* editor, where the document lives on one machine and
//! several people look at it, or where the browser cannot run the module. It is the one case the
//! local-first design does not serve, and it is deliberately the smaller half of the two: the
//! wasm binding needs nothing running, and this needs somebody to run it.
//!
//! **One document, named at startup, shared by every client.** A request sees the effect of the
//! request before it. What that does *not* solve is named rather than implied: there is no notion
//! of who edited what, no conflict resolution beyond the order requests arrive in, and no
//! authentication. It binds `127.0.0.1` and refuses every origin it was not told to allow, which
//! is the posture a single-user local tool can defend - anything else is a deployment, and a
//! deployment needs the things listed here that do not exist.
//!
//! # The protocol, which is deliberately small
//!
//! ```text
//! POST /rpc   {"command": {…}|null, "run": bool?, "order": "insertion"|"topological"?}
//!          →  200  the envelope
//!          →  400  {"error": "…"}      the *request* or the *call* was refused
//!
//! POST /rpc   {"catalogue": bool}
//!          →  200  the catalogue: one form per unit op, and the tools on request
//!          →  400  {"error": "…"}
//!
//! POST /mcp   one JSON-RPC message, in the current revision of MCP
//!          →  200  the message's answer
//!          →  202  accepted, and a notification has no answer
//!          →  400, 403, 404, 405, 406   see `crate::mcp_http`
//! ```
//!
//! A `command` the checker then refuses is **200 with `ok: false`** — the same three-way split the
//! MCP transport makes, in HTTP's own terms: a call that cannot take effect is a client error to
//! fix, while a document that cannot run is a state to show. `null` for `command` asks for the
//! envelope as it stands, which is what a client polls after somebody else's edit.
//!
//! **The catalogue is the one request that is not a call.** A front-end draws a palette and a form
//! per unit op, and neither is on the envelope, so a hosted editor has nowhere else to ask — but it
//! is a *read* of what this process loaded at startup and not an opening: nothing is handed a
//! document, nothing is replaced, and the session stays the one document `--flowsheet` named. It is
//! therefore its own body shape rather than a key on a call, because a body that asked for a
//! catalogue and an edit at once would have two answers and no way to choose between them.
//!
//! **`/mcp` is the same session behind a second door**, and the rules it adds are the ones its
//! protocol owns rather than ones about a flowsheet: a required protocol-version header that has to
//! agree with the body, a handshake this revision removed, and a status for each refusal. It is in
//! its own module because it is a different protocol speaking over the same socket rather than
//! another route of this one.
//!
//! **Hand-rolled, and the cost is the same one the MCP transport named.** HTTP/1.1 is a large
//! standard and this is a request line, a header block and a `Content-Length` body, answered with
//! `Connection: close`. No keep-alive, no chunked encoding, no compression, no HTTP/2 — a client
//! that wants any of those is a client this is not for. Requests are handled **one at a time** on
//! the accepting thread, so there is no shared mutable state and no lock; a slow client delays the
//! next one, which for one document and a handful of people is a trade rather than a bug.

use std::io::{BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::path::Path;

use azoth_core::{AzothError, Result};
use serde_json::{Value, json};

use crate::session::Session;

/// The request path the editor's calls go to. One route, because there is one kind of call.
const ROUTE: &str = "/rpc";

/// The request path an agent's calls go to — the same session, the same envelope, a second door.
///
/// **The two routes are one document.** An edit through `/mcp` is visible to `/rpc` on the next
/// request and the other way round, because there is one `Session` behind both; that is what makes
/// this a second door rather than a second server, and it is a test rather than a claim.
const MCP_ROUTE: &str = "/mcp";

/// The largest body accepted, in bytes.
///
/// A command is a few hundred bytes and an envelope request is two; a megabyte is generous and
/// still refuses a client that would otherwise make this process allocate without limit.
const MAX_BODY: usize = 1024 * 1024;

/// Serve the session until the process is killed.
///
/// **The bound address is printed**, because `--port 0` asks the operating system for a free port
/// and a caller has to be able to read which one it got - the test beside this does exactly that,
/// and a person starting a server wants the URL rather than a guess.
///
/// # Errors
/// A palette that does not load, a document that does not read, or a port that cannot be bound.
pub fn serve(
    flowsheet: &Path,
    palette_dir: &Path,
    port: u16,
    allow_origins: &[String],
) -> Result<()> {
    let mut session = Session::open(flowsheet, palette_dir, true)?;

    let listener = TcpListener::bind(("127.0.0.1", port))
        .map_err(|error| AzothError::invalid_input("port", format!("127.0.0.1:{port}: {error}")))?;
    let address = listener
        .local_addr()
        .map_err(|error| AzothError::invalid_input("port", error.to_string()))?;
    // The line a caller reads to learn the port, and the one a person reads to learn the URL.
    println!("azoth: serving http://{address}{ROUTE}");
    if allow_origins.is_empty() {
        println!("azoth: no origin is allowed, so only a non-browser client can call this");
    } else {
        println!("azoth: allowing origins {}", allow_origins.join(", "));
    }

    for stream in listener.incoming() {
        let stream =
            stream.map_err(|error| AzothError::invalid_input("socket", error.to_string()))?;
        if let Err(error) = handle(stream, &mut session, allow_origins) {
            // A connection that fails is a client's problem and not this process's: it is written
            // to stderr and the next one is accepted.
            eprintln!("azoth: {error}");
        }
    }
    Ok(())
}

/// One request, one response, and the connection closed.
fn handle(
    stream: TcpStream,
    session: &mut Session,
    allow_origins: &[String],
) -> std::result::Result<(), String> {
    let mut reader = BufReader::new(stream.try_clone().map_err(|e| e.to_string())?);
    let mut writer = stream;

    let request = match read_request(&mut reader) {
        Ok(request) => request,
        Err(error) => return write_json(&mut writer, 400, &json!({ "error": error }), None),
    };

    // The one header that decides whether a page in a browser is allowed to read the answer.
    let origin = request.header("origin");
    let allowed = origin
        .as_deref()
        .filter(|origin| allow_origins.iter().any(|allowed| allowed == origin));
    let cors = allowed.map(|origin| format!("Access-Control-Allow-Origin: {origin}\r\n"));

    if request.method == "OPTIONS" {
        // A preflight, which a browser sends before a `POST` that is not "simple". Its answer is
        // the policy and nothing else: an origin this server was not told to allow gets a refusal
        // *without* the header, which is what the browser needs to see to block the call.
        return match cors {
            Some(cors) => write_raw(&mut writer, 204, "", Some(&cors)),
            None => write_json(
                &mut writer,
                403,
                &json!({ "error": "origin not allowed" }),
                None,
            ),
        };
    }
    if request.path == MCP_ROUTE {
        // The agent's door, whose rules are its own because its protocol is: the status a refusal
        // is served with is the specification's, and the header checks are HTTP's rather than the
        // command model's. What is *not* here is anything about a flowsheet.
        return match crate::mcp_http::answer(&request, cors.as_deref(), session) {
            crate::mcp_http::Answer::Reply { status, body } => {
                write_json(&mut writer, status, &body, cors.as_deref())
            }
            // A notification: accepted, and there is nothing to answer with.
            crate::mcp_http::Answer::Accepted => write_raw(&mut writer, 202, "", cors.as_deref()),
        };
    }
    if request.method != "POST" || request.path != ROUTE {
        return write_json(
            &mut writer,
            404,
            &json!({ "error": format!("no route; every call is a POST to {ROUTE} or {MCP_ROUTE}") }),
            cors.as_deref(),
        );
    }

    let body: Value = match serde_json::from_slice(&request.body) {
        Ok(body) => body,
        Err(error) => {
            return write_json(
                &mut writer,
                400,
                &json!({ "error": format!("the request body is not JSON: {error}") }),
                cors.as_deref(),
            );
        }
    };

    let answer = dispatch(session, &body);
    match answer {
        Ok(envelope) => write_json(&mut writer, 200, &envelope, cors.as_deref()),
        Err(sentence) => write_json(
            &mut writer,
            400,
            &json!({ "error": sentence }),
            cors.as_deref(),
        ),
    }
}

/// One request body, as the call it names.
///
/// The body is `{"command": …, "run": …, "order": …}` and **every part of it is optional**, so the
/// three things a client can want are one shape: an edit, a run, and a look at the envelope as it
/// stands. `order` is applied first, because it is a property of the *next* run and a request that
/// set it and ran in the same call means the run it sets it for.
///
/// `{"catalogue": …}` is the fourth and the only one that is not a call, so it is read first and
/// alone: nothing else in the body means anything to an answer that is not an envelope.
fn dispatch(session: &mut Session, body: &Value) -> std::result::Result<Value, String> {
    let object = body
        .as_object()
        .ok_or_else(|| "the request body is not an object".to_string())?;

    if let Some(with_tools) = object.get("catalogue") {
        if object.len() != 1 {
            return Err("a `catalogue` request carries nothing else".to_string());
        }
        let with_tools = with_tools
            .as_bool()
            .ok_or_else(|| "`catalogue` takes true or false".to_string())?;
        return session.catalogue(with_tools);
    }

    if let Some(order) = object.get("order") {
        let name = order
            .as_str()
            .ok_or_else(|| "`order` takes a name".to_string())?;
        session.set_order(name)?;
    }

    let run = match object.get("run") {
        None => None,
        Some(Value::Bool(run)) => Some(*run),
        Some(_) => return Err("`run` takes true or false".to_string()),
    };

    match object.get("command") {
        // An edit, which runs the document as the invocation says unless this request says
        // otherwise.
        Some(command) if !command.is_null() => {
            let name = command
                .get("command")
                .and_then(Value::as_str)
                .ok_or_else(|| "a command names itself, as `{\"command\": \"…\"}`".to_string())?;
            session.call(name, command, run)
        }
        // No command: the envelope as it stands, which is what a client polls after somebody
        // else's edit.
        _ => session.settle(run),
    }
}

/// The parts of a request this server reads, and nothing else.
pub(crate) struct Request {
    pub(crate) method: String,
    pub(crate) path: String,
    headers: Vec<(String, String)>,
    pub(crate) body: Vec<u8>,
}

impl Request {
    /// A header's value, by name and **case-insensitively**: HTTP field names are, and a client
    /// that spells `mcp-method` is a conforming one.
    pub(crate) fn header(&self, name: &str) -> Option<String> {
        self.headers
            .iter()
            .find(|(key, _)| key.eq_ignore_ascii_case(name))
            .map(|(_, value)| value.clone())
    }
}

/// A request line, its headers and its body.
///
/// # Errors
/// A sentence for anything the parser will not guess at: a request line that is not three parts, a
/// header without a colon, a `Content-Length` that is not a number, a body longer than
/// [`MAX_BODY`], or a body that ends early.
pub(crate) fn read_request(reader: &mut impl BufRead) -> std::result::Result<Request, String> {
    let mut line = String::new();
    reader
        .read_line(&mut line)
        .map_err(|error| error.to_string())?;
    let mut parts = line.trim_end().split(' ');
    let (Some(method), Some(path), Some(version)) = (parts.next(), parts.next(), parts.next())
    else {
        return Err(format!("`{}` is not a request line", line.trim_end()));
    };
    if !version.starts_with("HTTP/") {
        return Err(format!("`{version}` is not an HTTP version"));
    }

    let mut headers = Vec::new();
    loop {
        let mut line = String::new();
        if reader
            .read_line(&mut line)
            .map_err(|error| error.to_string())?
            == 0
        {
            return Err("the headers ended without a blank line".to_string());
        }
        let line = line.trim_end();
        if line.is_empty() {
            break;
        }
        let (key, value) = line
            .split_once(':')
            .ok_or_else(|| format!("`{line}` is not a header"))?;
        headers.push((key.trim().to_string(), value.trim().to_string()));
    }

    let length: usize = headers
        .iter()
        .find(|(key, _)| key.eq_ignore_ascii_case("content-length"))
        .map_or(Ok(0), |(_, value)| {
            value
                .parse()
                .map_err(|_| format!("`{value}` is not a Content-Length"))
        })?;
    if length > MAX_BODY {
        return Err(format!(
            "a body of {length} bytes is over the {MAX_BODY} this accepts"
        ));
    }
    let mut body = vec![0_u8; length];
    reader
        .read_exact(&mut body)
        .map_err(|error| format!("the body ended early: {error}"))?;

    Ok(Request {
        method: method.to_string(),
        path: path.to_string(),
        headers,
        body,
    })
}

/// A JSON response, with the CORS header where one is allowed.
fn write_json(
    writer: &mut impl Write,
    status: u16,
    body: &Value,
    cors: Option<&str>,
) -> std::result::Result<(), String> {
    let text = serde_json::to_string(body).map_err(|error| error.to_string())?;
    write_raw(writer, status, &text, cors)
}

/// A response, with the headers this server always sends.
///
/// `Connection: close` on every one of them, which is what makes the missing keep-alive a decision
/// rather than a bug: the client is told, so it opens another connection.
fn write_raw(
    writer: &mut impl Write,
    status: u16,
    body: &str,
    cors: Option<&str>,
) -> std::result::Result<(), String> {
    let reason = match status {
        200 => "OK",
        202 => "Accepted",
        204 => "No Content",
        400 => "Bad Request",
        // The one status here no route of this server returns yet: `crate::mcp_http`'s payment
        // gate is dormant and answers nothing, and this is the status its seam would carry.
        402 => "Payment Required",
        403 => "Forbidden",
        404 => "Not Found",
        405 => "Method Not Allowed",
        406 => "Not Acceptable",
        _ => "Not Found",
    };
    let cors = cors.unwrap_or_default();
    // A response with nothing in it has no content to describe: a `Content-Type` on an empty body
    // is a claim about bytes that are not there.
    let content_type = if body.is_empty() {
        String::new()
    } else {
        "Content-Type: application/json\r\n".to_string()
    };
    let response = format!(
        "HTTP/1.1 {status} {reason}\r\n\
         {cors}\
         {content_type}\
         Content-Length: {}\r\n\
         Connection: close\r\n\
         \r\n\
         {body}",
        body.len()
    );
    writer
        .write_all(response.as_bytes())
        .and_then(|()| writer.flush())
        .map_err(|error| error.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn a_request() -> Request {
        Request {
            method: "POST".to_string(),
            path: MCP_ROUTE.to_string(),
            headers: Vec::new(),
            body: Vec::new(),
        }
    }

    /// **The payment gate refuses nothing, and that is measured rather than assumed.** The seam is
    /// dormant because this server does not charge — there is no token, no issuer and no quota — and
    /// a dormant branch with no test is the shape this repository's gates exist to catch.
    #[test]
    fn the_payment_gate_refuses_nothing_today() {
        assert!(crate::mcp_http::payment_refusal(&a_request()).is_none());
    }

    /// **And the status it would answer with is renderable**, which is the other half: the reason
    /// map is what decides whether a client is told "payment required" or "not found", and no route
    /// returns it yet.
    #[test]
    fn the_payment_required_status_is_renderable() {
        let mut written = Vec::new();
        write_raw(&mut written, 402, r#"{"error":"payment required"}"#, None).expect("it writes");
        let text = String::from_utf8(written).expect("it is text");
        assert!(
            text.starts_with("HTTP/1.1 402 Payment Required\r\n"),
            "{text}"
        );
        assert!(text.contains("Connection: close"), "{text}");
    }
}
