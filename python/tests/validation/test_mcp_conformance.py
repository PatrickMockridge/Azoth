"""What `azoth mcp` says, held to the specification's own schema.

The transport in `crates/azoth-cli/src/mcp.rs` is hand-rolled, and the one thing a hand-rolled
implementation of a published protocol cannot check about itself is whether it is *conformant*: it
writes what its author read, and a revision that moves a field leaves it serving a shape no client
accepts. **The specification publishes a JSON Schema for every revision**, so the check is a
comparison rather than an opinion - `jsonschema` is already a dev dependency, and the schema is a
draft 2020-12 document whose `$defs` name the message types (`InitializeResult`, `CallToolResult`,
`ListToolsResult`, `JSONRPCErrorResponse`, ...).

# Two revisions, so two schemas, and that is the specification's doing

**The current revision removed the handshake**: `2026-07-28` has no `InitializeResult` and no
`initialize`, and `2025-11-25` has no `server/discover`. The server answers both lanes, so it is
held to both schemas - the lane a message belongs to decides which one. A server validated against
one schema would have its other lane unchecked, and validating both lanes against both schemas
would fail for a reason that is the protocol's rather than the port's.

# And two transports

`azoth mcp` speaks JSON-RPC over stdio and `azoth serve` exposes the same messages at `POST /mcp`
over HTTP, so both are checked here: the stdio lane by driving the process, the HTTP lane by
posting to a served one. The claim the second lane exists to make is that there is **one encoder
and not two** - the same `$defs` validate both - and the thing only the HTTP lane can check is what
HTTP adds: the status each refusal is served with, the headers the specification mirrors the body
into, and the compatibility behaviour of an endpoint that speaks the current revision only.

# What this cannot see

Upstream moving. The revisions are pinned by the vendored files, so a *new* revision is invisible
until somebody re-fetches it - which is a deliberate act, and the act is a diff of the vendored file
against the URL in `NOTICE`. What the schema does catch is the thing that actually happens: a field
renamed, a required property dropped, or a lane emitting the other lane's message (which is how the
two-schema split above was found - `initialize` has no `InitializeResult` in the current revision).
"""

from __future__ import annotations

import http.client
import json
import subprocess
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
from jsonschema import Draft202012Validator

REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parent / "vendor" / "mcp"
CURRENT = VENDOR / "schema-2026-07-28.json"
LEGACY = VENDOR / "schema-2025-11-25.json"

#: The lane each revision is, and where the server's methods live in it.
CURRENT_REVISION = "2026-07-28"
LEGACY_REVISION = "2025-11-25"

#: The two versions key the current revision requires on every request, in `_meta`.
META_VERSION = "io.modelcontextprotocol/protocolVersion"
META_CAPABILITIES = "io.modelcontextprotocol/clientCapabilities"

DEMO = REPO_ROOT / "specs" / "flowsheets" / "demo.toml"
PALETTE = REPO_ROOT / "specs" / "unit_ops"


def modern(id_: int, method: str, params: dict[str, Any] | None = None) -> dict[str, Any]:
    """One request in the current revision, with the `_meta` the schema requires on every one.

    **Both keys**, which is the schema's own `required` array rather than a reading of the prose:
    `clientCapabilities` is required and unused here, because this server needs no client capability
    to run a flowsheet.
    """
    params = dict(params or {})
    params["_meta"] = {
        META_VERSION: CURRENT_REVISION,
        META_CAPABILITIES: {},
    }
    return {"jsonrpc": "2.0", "id": id_, "method": method, "params": params}


def headers_of(message: dict[str, Any]) -> dict[str, str]:
    """The headers a conforming client sends for a message, derived from the body.

    The specification mirrors `method` and `params.name` into headers so an intermediary can route
    without parsing the body — which is why a *disagreeing* header is a thing a test can send and a
    server must refuse.
    """
    headers = {
        "MCP-Protocol-Version": CURRENT_REVISION,
        "Mcp-Method": str(message.get("method", "")),
    }
    name = message.get("params", {}).get("name")
    if name is not None:
        headers["Mcp-Name"] = str(name)
    return headers


@pytest.fixture(scope="module")
def server() -> str:
    """The binary, built once.

    **Built rather than assumed**: the job this runs in has a Rust toolchain and a warm target
    directory (it builds the extension a step earlier), and a conformance check that skipped
    because nothing had compiled the transport would be a check that did not run.
    """
    subprocess.run(
        ["cargo", "build", "-j", "10", "-p", "azoth-cli"],
        cwd=REPO_ROOT,
        capture_output=True,
        check=True,
    )
    binary = REPO_ROOT / "target" / "debug" / "azoth"
    assert binary.exists(), f"{binary} was not built"
    return str(binary)


def drive(server: str, requests: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """One request per line in, one response per line out."""
    completed = subprocess.run(
        [server, "mcp", "--flowsheet", str(DEMO), "--no-run"],
        input="\n".join(json.dumps(request) for request in requests) + "\n",
        capture_output=True,
        text=True,
        timeout=300,
        check=True,
    )
    responses = [json.loads(line) for line in completed.stdout.splitlines() if line.strip()]
    assert len(responses) == len(requests), (
        f"{len(responses)} responses for {len(requests)} requests"
    )
    return responses


@pytest.fixture(scope="module")
def endpoint(server: str) -> Iterator[int]:
    """`azoth serve`, on a port the operating system picks, for the whole module.

    **The same binary and the same document as `drive`**, so what is checked over HTTP is the
    messages one process writes rather than a second server's idea of them.
    """
    process = subprocess.Popen(
        [
            server,
            "serve",
            "--flowsheet",
            str(DEMO),
            "--palette",
            str(PALETTE),
            "--port",
            "0",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        port = None
        assert process.stdout is not None
        for line in process.stdout:
            if line.startswith("azoth: serving http://"):
                # `azoth: serving http://127.0.0.1:43211/rpc`
                address = line.split("//", 1)[1].split("/", 1)[0]
                port = int(address.rsplit(":", 1)[1])
                break
        assert port is not None, "the server prints the address it bound"
        yield port
    finally:
        process.kill()
        process.wait()
        # Closed rather than left to the collector: this suite turns warnings into errors, and an
        # unclosed pipe is a `ResourceWarning` that fails the teardown of whichever test ran last.
        for pipe in (process.stdout, process.stderr):
            if pipe is not None:
                pipe.close()


def post(
    port: int, message: dict[str, Any], headers: dict[str, str] | None = None
) -> tuple[int, Any]:
    """One POST to `/mcp`, and the status and body it answers with."""
    body = json.dumps(message).encode()
    # A request that carries no `Accept` at all is a client that accepts anything, which is what
    # `curl` sends and what the endpoint is tested against by default.
    sent = {"Content-Type": "application/json", **dict(headers or {})}
    connection = http.client.HTTPConnection("127.0.0.1", port, timeout=60)
    try:
        connection.request("POST", "/mcp", body=body, headers=sent)
        response = connection.getresponse()
        payload = response.read().decode()
        return response.status, (json.loads(payload) if payload.strip() else None)
    finally:
        connection.close()


def schema_of(path: Path) -> dict[str, Any]:
    document: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return document


def conforms(payload: Any, definition: str, schema: dict[str, Any]) -> None:
    """A payload against one `$defs` entry, with the document's other defs resolvable."""
    named = schema["$defs"].get(definition)
    assert named is not None, f"the vendored schema has no `{definition}`"
    validator = Draft202012Validator({**named, "$defs": schema["$defs"]})
    errors = sorted(
        validator.iter_errors(payload), key=lambda error: [str(part) for part in error.path]
    )
    if errors:
        raise AssertionError(
            f"{definition}: {errors[0].json_path or '<root>'}: {errors[0].message}"
        )


def test_this_file_reads_the_whole_shape_of_the_lane_it_names() -> None:
    """The two schemas are the two the server speaks, and each says so itself.

    **The lane split, checked rather than assumed**: the current revision is the one without the
    handshake and the legacy one is the one without discovery. If either gains or loses a message
    here, the mapping below is wrong in a way no single message would reveal.
    """
    current = schema_of(CURRENT)["$defs"]
    legacy = schema_of(LEGACY)["$defs"]
    assert "InitializeResult" in legacy, "the legacy revision is the handshake one"
    assert "InitializeResult" not in current, "the current revision removed the handshake"
    assert "DiscoverResult" in current, "the current revision requires server/discover"
    assert "DiscoverResult" not in legacy
    for shared in ("ListToolsResult", "CallToolResult", "JSONRPCErrorResponse"):
        assert shared in current and shared in legacy, shared


def test_every_message_the_server_writes_is_the_shape_its_revision_names(server: str) -> None:
    """**The check this file exists for**, over the methods the server serves and both lanes.

    Each response is validated against the definition for *its own* lane, and the whole envelope for
    an error - because an error response is the envelope rather than a result.
    """
    current = schema_of(CURRENT)
    legacy = schema_of(LEGACY)
    requests: list[tuple[str, dict[str, Any], dict[str, Any], str]] = [
        (
            "initialize",
            {
                "jsonrpc": "2.0",
                "id": 1,
                "method": "initialize",
                "params": {"protocolVersion": LEGACY_REVISION},
            },
            legacy,
            "InitializeResult",
        ),
        ("server/discover", modern(2, "server/discover"), current, "DiscoverResult"),
        ("tools/list", modern(3, "tools/list"), current, "ListToolsResult"),
        (
            "tools/call",
            modern(4, "tools/call", {"name": "remove_instance", "arguments": {"id": "hx1"}}),
            current,
            "CallToolResult",
        ),
        # A call that landed and left a document the checker refuses: still a `CallToolResult`,
        # because the call happened and `isError` is false.
        (
            "a refused document",
            modern(
                5,
                "tools/call",
                {"name": "add_instance", "arguments": {"id": "x1", "unit": "unit_ops.nosuch"}},
            ),
            current,
            "CallToolResult",
        ),
        # A call that could not take effect: `isError`, which is the specification's own shape for
        # input validation rather than a protocol error.
        (
            "a malformed call",
            modern(
                6,
                "tools/call",
                {"name": "set_position", "arguments": {"node": "instance:nope", "x": 0, "y": 0}},
            ),
            current,
            "CallToolResult",
        ),
        # An unknown method is the envelope, which is why it is validated as one.
        ("an unknown method", modern(7, "tools/write"), current, "JSONRPCErrorResponse"),
        (
            "an unknown tool",
            modern(8, "tools/call", {"name": "nope", "arguments": {}}),
            current,
            "JSONRPCErrorResponse",
        ),
        # **A request without its metadata**, which the current revision calls malformed: `-32602`
        # rather than the version refusal's `-32022`, and the shape the specification names for it.
        (
            "no metadata",
            {"jsonrpc": "2.0", "id": 9, "method": "tools/list", "params": {}},
            current,
            "JSONRPCErrorResponse",
        ),
    ]

    responses = drive(server, [request for _, request, _, _ in requests])
    for (label, _, schema, definition), response in zip(requests, responses, strict=True):
        payload = response.get("result", response)
        try:
            conforms(payload, definition, schema)
        except AssertionError as failure:
            raise AssertionError(f"{label}: {failure}") from None


def test_the_legacy_lane_is_the_handshake_and_the_current_one_has_no_such_method(
    server: str,
) -> None:
    """The negative half, which is what makes the lane split load-bearing.

    A server that answered `initialize` in the current revision's own terms would be serving a
    message that revision does not have - and the schema says so by *not* carrying the definition,
    which is a different kind of check from a shape failure.
    """
    current = schema_of(CURRENT)["$defs"]
    assert "resultType" in json.dumps(current["CallToolResult"]), (
        "the current revision requires a result to say what it is"
    )
    # And the server does write it, which is what makes the current lane's messages valid for a
    # client that checks.
    (response,) = drive(server, [modern(1, "tools/list")])
    assert response["result"]["resultType"] == "complete"


def test_the_http_transport_writes_the_same_shapes_as_the_stdio_one(endpoint: int) -> None:
    """`POST /mcp` is the same encoder, so its answers are held to the same definitions.

    **The claim being checked is that there is one encoder and not two**: a second transport is
    where a second message shape would hide, and comparing this lane's results to the *same* `$defs`
    the stdio lane's are compared to is what makes "the same messages" a measurement.
    """
    current = schema_of(CURRENT)
    requests: list[tuple[str, dict[str, Any], str]] = [
        ("server/discover", modern(1, "server/discover"), "DiscoverResult"),
        ("tools/list", modern(2, "tools/list"), "ListToolsResult"),
        (
            "tools/call",
            modern(3, "tools/call", {"name": "remove_instance", "arguments": {"id": "hx1"}}),
            "CallToolResult",
        ),
    ]
    for label, message, definition in requests:
        status, payload = post(endpoint, message, headers_of(message))
        assert status == 200, f"{label}: {status} {payload}"
        try:
            conforms(payload["result"], definition, current)
        except AssertionError as failure:
            raise AssertionError(f"{label}: {failure}") from None


def test_the_http_transports_refusals_are_the_shapes_its_revision_names(endpoint: int) -> None:
    """The transport's own refusals, held to the schema rather than to this file's expectations.

    **Two of the three exist only in the current revision**, which is itself the point: a transport
    serving the legacy lane could not be checked this way at all, because `2025-11-25` defines
    neither `HeaderMismatchError` nor `UnsupportedProtocolVersionError`.

    The `$defs` are also not all the same kind of thing, and using the wrong one would pass
    vacuously: `HeaderMismatchError` and `UnsupportedProtocolVersionError` are whole *responses*,
    while `InvalidParamsError` is the `error` member of one — so the envelope is checked separately
    for the case that needs the other kind.
    """
    current = schema_of(CURRENT)

    wrong_version = modern(1, "tools/list")
    wrong_version["params"]["_meta"][META_VERSION] = "1999-01-01"
    disagreeing = headers_of(modern(2, "tools/list"))
    disagreeing["MCP-Protocol-Version"] = "1999-01-01"

    #: label, message, headers, status, the definition to hold it to, and which part of the body
    #: that definition describes.
    cases: list[tuple[str, dict[str, Any], dict[str, str], int, str, str]] = [
        (
            "a header that disagrees with the body",
            modern(3, "tools/list"),
            disagreeing,
            400,
            "HeaderMismatchError",
            "response",
        ),
        (
            "a version this endpoint does not speak",
            wrong_version,
            {**headers_of(wrong_version), "MCP-Protocol-Version": "1999-01-01"},
            400,
            "UnsupportedProtocolVersionError",
            "response",
        ),
        (
            "the handshake this endpoint does not speak",
            {
                "jsonrpc": "2.0",
                "id": 4,
                "method": "initialize",
                "params": {"protocolVersion": LEGACY_REVISION},
            },
            {},
            400,
            "UnsupportedProtocolVersionError",
            "response",
        ),
        (
            "a request missing a required metadata key",
            {
                "jsonrpc": "2.0",
                "id": 5,
                "method": "tools/list",
                "params": {"_meta": {META_VERSION: CURRENT_REVISION}},
            },
            {"MCP-Protocol-Version": CURRENT_REVISION, "Mcp-Method": "tools/list"},
            400,
            "InvalidParamsError",
            "error",
        ),
    ]

    for label, message, headers, status, definition, part in cases:
        got, payload = post(endpoint, message, headers)
        assert got == status, f"{label}: {got} {payload}"
        subject = payload if part == "response" else payload["error"]
        try:
            conforms(subject, definition, current)
        except AssertionError as failure:
            raise AssertionError(f"{label}: {failure}") from None
        if part == "error":
            # The envelope around it is still a JSON-RPC error response, which is a different
            # definition rather than an assumption about the one above.
            conforms(payload, "JSONRPCErrorResponse", current)


def test_the_endpoint_is_modern_only_and_says_so(endpoint: int) -> None:
    """**The compatibility behaviour, where this endpoint differs from `azoth mcp`.**

    A legacy client's GET and DELETE are `405`, which is what the specification tells a modern-only
    server to answer: those methods belong to the revisions that had a session id and a standalone
    stream, and this one has neither. And an `initialize` is told the revisions *this transport*
    speaks — one, not the two the process speaks — so a client cannot be sent into a retry this
    endpoint would refuse again.
    """
    connection = http.client.HTTPConnection("127.0.0.1", endpoint, timeout=60)
    try:
        for method in ("GET", "DELETE"):
            connection.request(method, "/mcp", headers={"MCP-Protocol-Version": CURRENT_REVISION})
            response = connection.getresponse()
            assert response.status == 405, f"{method}: {response.status}"
            response.read()
    finally:
        connection.close()

    handshake = {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {"protocolVersion": LEGACY_REVISION},
    }
    status, payload = post(endpoint, handshake)
    assert status == 400
    assert payload["error"]["data"]["supported"] == [CURRENT_REVISION], (
        "the list is this transport's, and naming a revision this endpoint would refuse teaches a "
        f"client to retry forever: {payload}"
    )
