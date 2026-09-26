//! `POST /mcp` — the tool schema over Streamable HTTP.
//!
//! **The third door and not a third surface.** The messages are [`crate::mcp`]'s own — the same
//! `Protocol`, the same tools, the same envelope — so an agent driving this and an agent driving
//! `azoth mcp` over stdio get the same answers to the same calls. What this file adds is everything
//! HTTP owns and stdio does not: a status for each refusal, the headers the specification mirrors
//! the body into, and the two things a *modern-only* transport owes a client that speaks an older
//! revision.
//!
//! **The endpoint is on `azoth serve`, so the document is shared.** An agent over HTTP and the
//! browser editor edit one `Workspace`, which is the whole reason this is a route rather than a
//! second binary: two transports over one document is the middleware's claim, and two processes
//! would be two documents.
//!
//! # Modern only, and what that decides
//!
//! This endpoint speaks `2026-07-28`. It is not a shortcut but a fit: that revision removed
//! protocol-level sessions and made every request carry its own metadata, so it holds nothing
//! between requests — and neither does this, which is why a `Protocol` is made per request below.
//! The legacy lane's HTTP shape is a session id, a standalone GET stream and resumable streams, all
//! of which exist to carry server-initiated messages, and this server sends none.
//!
//! So an older client's traffic is answered as the specification tells a modern-only server to
//! answer it: `405` for the GET and DELETE of the revisions that had them, `Mcp-Session-Id` and
//! `Last-Event-ID` ignored rather than honoured, and an `initialize` that names the revisions this
//! transport *does* speak. **`azoth mcp` over stdio still speaks both**, so nothing is lost — the
//! handshake is one process away.
//!
//! # The order of the checks, and why it is this order
//!
//! Origin, then payment, then the method, then the body, then the headers, then the protocol. The
//! two security checks come first because they need only headers; the protocol comes last because
//! everything before it can refuse without parsing a document. The one the specification pins and
//! this file exists to get right is the header check: `MCP-Protocol-Version` must equal the body's
//! `_meta` version, `Mcp-Method` its `method`, and `Mcp-Name` its `params.name` — *"so that
//! intermediaries can route and inspect requests without parsing the body"*, with the body as the
//! source of truth and `-32020 HeaderMismatch` when the two disagree. A gateway routing on a header
//! while this process executes the body is exactly the disagreement that code is for.

use serde_json::{Value, json};

use crate::mcp::{self, Protocol};
use crate::serve::Request;
use crate::session::Session;

/// The protocol revisions **this transport** speaks.
///
/// One, where the process speaks two: `azoth mcp` answers the handshake over stdio and this does
/// not. **The list is what a legacy client is told**, so naming `2025-11-25` here would send it
/// into a retry this endpoint can never accept — the honest answer names what would work (the
/// current revision) and the message says where the other lane lives.
const SUPPORTED: [&str; 1] = ["2026-07-28"];

/// The request-metadata key the version header must agree with.
const META_VERSION: &str = "io.modelcontextprotocol/protocolVersion";

/// The headers the specification mirrors from the body, and what each one is checked against.
const HEADER_VERSION: &str = "MCP-Protocol-Version";
const HEADER_METHOD: &str = "Mcp-Method";
const HEADER_NAME: &str = "Mcp-Name";

/// The requests whose `Mcp-Name` is required, in the specification's own list.
///
/// The other two are not served by this server — but the *header* requirement is the transport's
/// and is answered before the method is judged, so a client sending `resources/read` without a
/// `Mcp-Name` is told what is missing rather than that the method does not exist.
const NAME_REQUIRED: [&str; 3] = ["tools/call", "resources/read", "prompts/get"];

/// JSON-RPC's `Parse error`, which is the one refusal with no id to answer to.
const PARSE_ERROR: i64 = -32700;

/// `HeaderMismatch`, the specification's own code for this transport's one new refusal.
const HEADER_MISMATCH: i64 = -32020;

/// `UnsupportedProtocolVersion`, for the handshake this revision removed.
const UNSUPPORTED_VERSION: i64 = -32022;

/// What the transport decided, for `serve.rs` to put on the wire.
///
/// **The status is decided here and written there**, so every rule about which status a refusal
/// gets is in one file and the socket layer stays a socket layer.
pub(crate) enum Answer {
    /// A status and a body.
    Reply { status: u16, body: Value },
    /// `202 Accepted`, with nothing in it: a notification has no reply, and answering one puts a
    /// message on the wire the client is not reading.
    Accepted,
}

/// One message, as the answer it gets.
pub(crate) fn answer(request: &Request, cors: Option<&str>, session: &mut Session) -> Answer {
    // The origin first, and this transport is stricter than `/rpc` on purpose. The specification:
    // *"Servers MUST validate the `Origin` header on all incoming connections to prevent DNS
    // rebinding attacks. If the `Origin` header is present and invalid, servers MUST respond with
    // HTTP 403 Forbidden."* `/rpc` withholds the CORS header instead — which is what a browser
    // reads as "no" — and is left exactly as it is, because the editor's door is defended by the
    // browser while this one may be reached by a client that is not one.
    if request.header("origin").is_some() && cors.is_none() {
        return self_refusal(403, "origin not allowed".to_string());
    }

    // The payment gate, which is dormant: see `payment_refusal`.
    if let Some(refusal) = payment_refusal(request) {
        return self_refusal(refusal.status, refusal.sentence);
    }

    // The method. GET, DELETE and the standalone SSE stream they opened were removed by the
    // revision this endpoint speaks, and the specification tells a server in this position to
    // answer them `405` rather than to pretend they are routes.
    if request.method != "POST" {
        return self_refusal(
            405,
            format!(
                "the MCP endpoint takes POST, and `{}` is what a revision this server does not \
                 speak uses",
                request.method
            ),
        );
    }

    // `Accept`, which a JSON-only server has exactly one use for: to refuse honestly. A client that
    // says it cannot read `application/json` has told us it cannot read the only answer this
    // endpoint writes, and sending it anyway would be a reply it is entitled to discard.
    if !accepts_json(request) {
        let accept = request.header("accept").unwrap_or_default();
        return self_refusal(
            406,
            format!("this endpoint answers `application/json`, and `{accept}` does not accept it"),
        );
    }

    let body: Value = match serde_json::from_slice(&request.body) {
        Ok(body) => body,
        Err(error) => {
            return Answer::Reply {
                status: 400,
                body: mcp::error_response(&Value::Null, PARSE_ERROR, &error.to_string(), None),
            };
        }
    };

    // **A notification is accepted and not answered**, and the header rules below are the
    // specification's for *requests* — it defines none for a notification POST, so checking them
    // here would refuse a message it says to accept. One predicate decides both, and the `None`
    // from `Protocol::handle` below is the same rule arriving at the same answer.
    if !mcp::is_notification(&body) {
        // **The handshake first, and this order is the whole of it.** A legacy client could not
        // have sent the headers below — its revision never defined them, and its version travelled
        // in the handshake's own `params` — so judging it by the header rule would answer a client
        // that did nothing wrong with a complaint about a header it had no way to send. The
        // compatibility answer comes first for that reason, and the header rule then applies to the
        // messages that are trying to be modern.
        if body.get("method").and_then(Value::as_str) == Some("initialize") {
            return Answer::Reply {
                status: 400,
                body: handshake_refused(&body),
            };
        }
        if let Err(sentence) = headers_agree(request, &body) {
            return Answer::Reply {
                status: 400,
                body: mcp::error_response(
                    &body.get("id").cloned().unwrap_or(Value::Null),
                    HEADER_MISMATCH,
                    &sentence,
                    None,
                ),
            };
        }
    }

    // **A fresh `Protocol` per request, because this revision is stateless**: nothing may be
    // inferred from a request that came before, and there is nothing to infer — the lane is the
    // metadata the request carries. What it holds for stdio is the lane a client opened in, and a
    // process is a scope while an endpoint is not.
    let mut protocol = Protocol::new();
    match protocol.handle(session, &body) {
        Some(response) => Answer::Reply {
            status: status_of(&response),
            body: response,
        },
        None => Answer::Accepted,
    }
}

/// The status a JSON-RPC answer is served with.
///
/// **The specification pins one of them and this is the only place the mapping lives**: a method
/// this server does not serve is `404 Not Found`, and every other refusal it names — a header that
/// disagrees with the body, a version it does not speak, a request missing its `_meta` — is `400
/// Bad Request`. A *result* is `200`, including the two states a caller has to read rather than
/// correct: a call that could not take effect (`isError`) and a document the checker refuses
/// (`ok: false`), because in both of them the call happened.
fn status_of(response: &Value) -> u16 {
    match response.pointer("/error/code").and_then(Value::as_i64) {
        Some(mcp::METHOD_NOT_FOUND) => 404,
        Some(_) => 400,
        None => 200,
    }
}

/// Whether the headers agree with the body, as the specification requires them to.
///
/// Three headers, three comparisons, and every refusal is the same code: `-32020 HeaderMismatch`,
/// which the schema describes as *"returned when a server rejects a request because the values in
/// the HTTP headers do not match the corresponding values in the request body, or because required
/// headers are missing or malformed"*. Both halves are here — a missing header and a disagreeing
/// one — and the sentences differ because a client fixes them differently.
fn headers_agree(request: &Request, body: &Value) -> std::result::Result<(), String> {
    let method = body
        .get("method")
        .and_then(Value::as_str)
        .unwrap_or_default();
    // Two steps rather than one pointer, because the key has a `/` in it and a JSON Pointer would
    // have to escape it — an escaping that is easy to get wrong and impossible to see.
    let version = body
        .pointer("/params/_meta")
        .and_then(|meta| meta.get(META_VERSION))
        .and_then(Value::as_str);

    agree(request, HEADER_VERSION, version)?;
    agree(request, HEADER_METHOD, Some(method))?;
    if NAME_REQUIRED.contains(&method) {
        let name = body
            .pointer("/params/name")
            .or_else(|| body.pointer("/params/uri"))
            .and_then(Value::as_str);
        agree(request, HEADER_NAME, name)?;
    }
    Ok(())
}

/// One header against the body's own value for it, with the value decoded first.
fn agree(request: &Request, header: &str, wanted: Option<&str>) -> std::result::Result<(), String> {
    let Some(given) = request.header(header) else {
        return Err(format!("`{header}` is required for this request"));
    };
    // **Decoded before it is compared**, which the specification requires of a server that inspects
    // these values: a header carrying Base64 is what a name outside the header-safe set travels as,
    // and comparing the encoded text against the body would refuse a conforming client calling a
    // tool this server does serve.
    let given = sentinel_decode(&given).map_err(|sentence| format!("`{header}` {sentence}"))?;
    let Some(wanted) = wanted else {
        return Err(format!(
            "`{header}` says `{given}` and the body names nothing to compare it against"
        ));
    };
    if given != wanted {
        return Err(format!(
            "`{header}` says `{given}` and the body says `{wanted}`"
        ));
    }
    Ok(())
}

/// The `400` an `initialize` gets here.
///
/// **The list is this transport's and the message says where the other lane is**, because a legacy
/// client has no fall-forward: a `2025-11-25` client told only `400` has nothing to do next, and
/// one told the versions would otherwise retry with a revision this endpoint will refuse again.
/// `data.requested` comes from the handshake's own `params.protocolVersion` — the schema requires
/// `requested` and `supported` together, so when a client names no version the list goes in the
/// message instead and `data` is left out rather than filled with a string that is not a version.
fn handshake_refused(body: &Value) -> Value {
    let asked = body
        .pointer("/params/protocolVersion")
        .and_then(Value::as_str);
    let message = match asked {
        Some(asked) => format!(
            "`initialize` asked for `{asked}`, which this endpoint does not speak: `{}` removed \
             the handshake and this door serves the current revision. `azoth mcp` speaks both \
             revisions over stdio",
            SUPPORTED[0]
        ),
        None => format!(
            "`initialize` is the handshake `{}` removed, and this door serves `{}` only; \
             `azoth mcp` speaks both revisions over stdio",
            SUPPORTED[0], SUPPORTED[0]
        ),
    };
    mcp::error_response(
        &body.get("id").cloned().unwrap_or(Value::Null),
        UNSUPPORTED_VERSION,
        &message,
        asked.map(|asked| json!({ "requested": asked, "supported": SUPPORTED })),
    )
}

/// A refusal this transport makes in its own terms rather than the protocol's.
///
/// The body is `azoth serve`'s own shape — `{"error": "…"}` — because these refusals are about the
/// request rather than about a message: an origin, a method, a content type. A JSON-RPC error body
/// would be this file borrowing the protocol's shape for something the protocol did not decide.
fn self_refusal(status: u16, sentence: String) -> Answer {
    Answer::Reply {
        status,
        body: json!({ "error": sentence }),
    }
}

/// A refusal decided before the protocol layer sees the request.
pub(crate) struct Refusal {
    pub(crate) status: u16,
    pub(crate) sentence: String,
}

/// **The payment gate, and it refuses nothing.**
///
/// This server does not charge: there is no token, no issuer and no quota, so there is no request
/// that could make a payment due and `None` is not a placeholder for a decision made elsewhere —
/// it is the whole of the decision.
///
/// It exists because the transport has to be able to answer *payment required* rather than *bad
/// request* the day this is deployed as a tokenized service, and the two are different instructions
/// to a client: `400` says the request was wrong and a retry will not help, while `402` says the
/// request was fine and there is an account to settle. The status is already in `serve.rs`'s reason
/// map, the body shape is [`self_refusal`]'s, and the call site is above — so what a tokenized
/// deployment adds is the credential in this function and nothing else.
///
/// **Dormant, and still measured**: `the_payment_gate_refuses_nothing_today` calls it, and
/// `a_refusal_is_the_servers_own_shape` renders a `402` through the same body a live one would
/// take, so the status and the reason phrase are on the wire in a test rather than being discovered
/// by a deployment.
pub(crate) fn payment_refusal(_request: &Request) -> Option<Refusal> {
    None
}

/// Whether a client can read the one content type this endpoint writes.
///
/// An absent header accepts anything, which is what makes `curl` work; a wildcard does too.
fn accepts_json(request: &Request) -> bool {
    let Some(accept) = request.header("accept") else {
        return true;
    };
    accept.contains("application/json")
        || accept.contains("application/*")
        || accept.contains("*/*")
}

/// The marker pair the specification puts around a Base64 header value.
const SENTINEL_PREFIX: &str = "=?base64?";
const SENTINEL_SUFFIX: &str = "?=";

/// A header's value, decoded if it was encoded.
///
/// A tool name is only *should*-constrained to header-safe characters, so a name with a non-ASCII
/// character, a control character or surrounding whitespace travels Base64-encoded between the
/// markers — and a plain value that happens to match the sentinel is encoded too, which is what
/// makes the markers unambiguous.
fn sentinel_decode(value: &str) -> std::result::Result<String, String> {
    let trimmed = value.trim();
    let Some(rest) = trimmed.strip_prefix(SENTINEL_PREFIX) else {
        return Ok(trimmed.to_string());
    };
    let Some(encoded) = rest.strip_suffix(SENTINEL_SUFFIX) else {
        return Err(format!(
            "starts the encoded-value marker (`{SENTINEL_PREFIX}`) and does not end it (`{SENTINEL_SUFFIX}`)"
        ));
    };
    base64_decode(encoded).ok_or_else(|| format!("carries `{encoded}`, which is not Base64"))
}

/// Standard Base64, with padding, and `None` for anything that is not.
///
/// **Hand-rolled because a dependency would be one crate for one comparison**: the alphabet is 64
/// characters and the padding is `=`, this file only ever *decodes*, and the specification's own
/// five encoding examples are what the unit test compares against — so the arithmetic is checked
/// against the text that defines it rather than against a second implementation of it.
fn base64_decode(text: &str) -> Option<String> {
    let mut bits = 0_u32;
    let mut held = 0_u32;
    let mut out = Vec::with_capacity(text.len() / 4 * 3);
    for (at, byte) in text.bytes().enumerate() {
        if byte == b'=' {
            // Padding ends the value and carries no bits, and nothing but padding may follow it.
            if text.bytes().skip(at).any(|byte| byte != b'=') {
                return None;
            }
            break;
        }
        let value = match byte {
            b'A'..=b'Z' => byte - b'A',
            b'a'..=b'z' => byte - b'a' + 26,
            b'0'..=b'9' => byte - b'0' + 52,
            b'+' => 62,
            b'/' => 63,
            _ => return None,
        } as u32;
        bits = (bits << 6) | value;
        held += 6;
        if held >= 8 {
            held -= 8;
            out.push((bits >> held) as u8);
            bits &= (1 << held) - 1;
        }
    }
    String::from_utf8(out).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The decoder against the text that defines it.** The specification publishes a table of
    /// five values and their encodings, so the arithmetic here is checked against the examples
    /// rather than against a second implementation of Base64 — and the two lines after it are the
    /// cases a table of encodings does not cover: a value that was never encoded, and one that
    /// starts the marker without ending it.
    #[test]
    fn the_sentinel_decodes_what_the_specification_says_it_does() {
        for (encoded, plain) in [
            ("us-west1", "us-west1"),
            ("=?base64?SGVsbG8sIOS4lueVjA==?=", "Hello, 世界"),
            ("=?base64?IHBhZGRlZCA=?=", " padded "),
            ("=?base64?bGluZTEKbGluZTI=?=", "line1\nline2"),
            ("=?base64?PT9iYXNlNjQ/bGl0ZXJhbD89?=", "=?base64?literal?="),
        ] {
            assert_eq!(
                sentinel_decode(encoded).expect("it decodes"),
                plain,
                "{encoded}"
            );
        }

        // Not encoded at all: a plain value is itself, which is what makes the comparison below it
        // a comparison of two names rather than of two encodings.
        assert_eq!(
            sentinel_decode("remove_instance").expect("it passes"),
            "remove_instance"
        );
        assert_eq!(sentinel_decode("  spaced  ").expect("it trims"), "spaced");

        // And the two ways an encoding can be malformed. Both are refusals rather than a silent
        // fallback to the raw text, because falling back would compare an encoding against a name.
        assert!(sentinel_decode("=?base64?cmVtb3ZlX2luc3RhbmNl").is_err());
        assert!(sentinel_decode("=?base64?not base64!?=").is_err());
    }

    /// **The dormant 402 is written in this server's own shape**, which is the half of the seam a
    /// live one would keep: the status is the instruction to the client and the body is the same
    /// `{"error": "…"}` every other transport refusal uses, so a tokenized deployment adds the
    /// credential and nothing else.
    #[test]
    fn a_payment_refusal_is_written_in_the_servers_own_shape() {
        let Answer::Reply { status, body } = self_refusal(402, "payment required".to_string())
        else {
            panic!("a refusal is a reply and not an acceptance");
        };
        assert_eq!(status, 402);
        assert_eq!(body, json!({ "error": "payment required" }));
    }
}
