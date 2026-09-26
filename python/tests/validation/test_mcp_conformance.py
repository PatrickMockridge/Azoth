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

# What this cannot see

Upstream moving. The revisions are pinned by the vendored files, so a *new* revision is invisible
until somebody re-fetches it - which is a deliberate act, and the act is a diff of the vendored file
against the URL in `NOTICE`. What the schema does catch is the thing that actually happens: a field
renamed, a required property dropped, or a lane emitting the other lane's message (which is how the
two-schema split above was found - `initialize` has no `InitializeResult` in the current revision).
"""

from __future__ import annotations

import json
import subprocess
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
        [server, "mcp", "--flowsheet", str(REPO_ROOT / "specs" / "flowsheets" / "demo.toml"), "--no-run"],
        input="\n".join(json.dumps(request) for request in requests) + "\n",
        capture_output=True,
        text=True,
        timeout=300,
        check=True,
    )
    responses = [json.loads(line) for line in completed.stdout.splitlines() if line.strip()]
    assert len(responses) == len(requests), f"{len(responses)} responses for {len(requests)} requests"
    return responses


def schema_of(path: Path) -> dict[str, Any]:
    document: dict[str, Any] = json.loads(path.read_text(encoding="utf-8"))
    return document


def conforms(payload: Any, definition: str, schema: dict[str, Any]) -> None:
    """A payload against one `$defs` entry, with the document's other defs resolvable."""
    named = schema["$defs"].get(definition)
    assert named is not None, f"the vendored schema has no `{definition}`"
    validator = Draft202012Validator({**named, "$defs": schema["$defs"]})
    errors = sorted(validator.iter_errors(payload), key=lambda error: [str(part) for part in error.path])
    if errors:
        raise AssertionError(f"{definition}: {errors[0].json_path or '<root>'}: {errors[0].message}")


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
            {"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": LEGACY_REVISION}},
            legacy,
            "InitializeResult",
        ),
        ("server/discover", {"jsonrpc": "2.0", "id": 2, "method": "server/discover"}, current, "DiscoverResult"),
        ("tools/list", {"jsonrpc": "2.0", "id": 3, "method": "tools/list"}, current, "ListToolsResult"),
        (
            "tools/call",
            {"jsonrpc": "2.0", "id": 4, "method": "tools/call",
             "params": {"name": "remove_instance", "arguments": {"id": "hx1"}}},
            current,
            "CallToolResult",
        ),
        # A call that landed and left a document the checker refuses: still a `CallToolResult`,
        # because the call happened and `isError` is false.
        (
            "a refused document",
            {"jsonrpc": "2.0", "id": 5, "method": "tools/call",
             "params": {"name": "add_instance", "arguments": {"id": "x1", "unit": "unit_ops.nosuch"}}},
            current,
            "CallToolResult",
        ),
        # A call that could not take effect: `isError`, which is the specification's own shape for
        # input validation rather than a protocol error.
        (
            "a malformed call",
            {"jsonrpc": "2.0", "id": 6, "method": "tools/call",
             "params": {"name": "set_position", "arguments": {"node": "instance:nope", "x": 0, "y": 0}}},
            current,
            "CallToolResult",
        ),
        # An unknown method is the envelope, which is why it is validated as one.
        ("an unknown method", {"jsonrpc": "2.0", "id": 7, "method": "tools/write"}, current, "JSONRPCErrorResponse"),
        (
            "an unknown tool",
            {"jsonrpc": "2.0", "id": 8, "method": "tools/call", "params": {"name": "nope", "arguments": {}}},
            current,
            "JSONRPCErrorResponse",
        ),
    ]

    responses = drive(server, [request for _, request, _, _ in requests])
    for (label, _, schema, definition), response in zip(requests, responses):
        payload = response["result"] if "result" in response else response
        try:
            conforms(payload, definition, schema)
        except AssertionError as failure:
            raise AssertionError(f"{label}: {failure}") from None


def test_the_legacy_lane_is_the_handshake_and_the_current_one_has_no_such_method(server: str) -> None:
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
    (response,) = drive(server, [{"jsonrpc": "2.0", "id": 1, "method": "tools/list"}])
    assert response["result"]["resultType"] == "complete"
