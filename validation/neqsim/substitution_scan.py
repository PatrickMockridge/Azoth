#!/usr/bin/env python3
"""Enumerate the substitution idiom of equinor/neqsim#4080.

    python3 scan.py <src-root>            # the definition below
    python3 scan.py <src-root> --wide     # also counts declaration initialisers

DEFINITION. One site is an assignment `name = <numeric literal>;` such that:

  1. it sits in the body of an `if` whose condition tests a NaN/Infinity *and* compares
     against zero, and the assigned `name` appears in that condition — so the guard is
     testing the very value it replaces;
  2. it is at depth 1 of that `if` (an assignment inside a nested block belongs to the
     nested guard, not this one, and is counted once);
  3. the `if` body contains no `throw` — a guard that refuses is not this idiom;
  4. comments are stripped first, so a trailing `// comment` does not disqualify a site.

`--wide` additionally counts `double name = <literal>;` declarations inside such a body.
Test sources are excluded.

MEASURED. Against NeqSim master `cf7e1c7ee9597e5df7e2b315cec4314b55c17351` (3.23.0), on
`src/main/java`, this returns **58 sites in 40 files**; `--wide` returns the same, because
no depth-1 declaration initialiser survives the nested-block rule. The capture is
`captures/substitution_scan.txt`.

    git clone --depth 1 https://github.com/equinor/neqsim
    python3 substitution_scan.py neqsim/src/main/java

The count is a *scope* statement, not the evidence: the list mixes legitimate named
defaults with substitutions that hide a failure, and separating those two is the
judgement the exemplar table in the report makes by hand. What this fixes is only that
the number is now reproducible from a definition stated here rather than depending on an
unrecorded one.
"""

import pathlib
import re
import sys

IF = re.compile(r"\bif\s*\(([^;{}]*)\)\s*\{", re.S)
ZERO = re.compile(r"[<>]=?\s*0(?:\.0*)?\b")
NAN = re.compile(r"isNaN|isInfinite|isFinite")
ASSIGN = re.compile(
    r"(?<![\w.])([A-Za-z_][A-Za-z0-9_]*(?:\[[^\]]*\])?(?:\.[A-Za-z_][A-Za-z0-9_]*)*)"
    r"\s*=\s*([-+]?\d[\d_.eE+\-]*[dDfFlL]?)\s*;"
)
DECL = re.compile(
    r"\b(?:double|float|int|long)\s+([A-Za-z_][A-Za-z0-9_]*)"
    r"\s*=\s*([-+]?\d[\d_.eE+\-]*[dDfFlL]?)\s*;"
)


def strip_comments(text):
    """Blank out comments and keep string literals intact, preserving offsets."""
    out, i, n = [], 0, len(text)
    while i < n:
        c = text[i]
        if c == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            out.append(text[i : j + 1])
            i = j + 1
        elif text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif text.startswith("/*", i):
            j = text.find("*/", i)
            j = n if j < 0 else j + 2
            out.append("".join(ch if ch == "\n" else " " for ch in text[i:j]))
            i = j
        else:
            out.append(c)
            i += 1
    return "".join(out)


def body_at_depth_one(text, start):
    """The depth-1 body of the block opening at `start`, nested blocks blanked out."""
    i, depth, out = start, 1, []
    while i < len(text) and depth:
        c = text[i]
        if c == "{":
            depth += 1
            if depth == 2:
                j, d = i + 1, 1
                while j < len(text) and d:
                    if text[j] == "{":
                        d += 1
                    elif text[j] == "}":
                        d -= 1
                    j += 1
                out.append(" " * (j - i))
                i = j
                depth = 1
                continue
        elif c == "}":
            depth -= 1
            if depth == 0:
                break
        out.append(c)
        i += 1
    return "".join(out)


def scan(path, wide):
    text = strip_comments(path.read_text(errors="replace"))
    hits = []
    for m in IF.finditer(text):
        cond = m.group(1)
        if not NAN.search(cond) or not ZERO.search(cond):
            continue
        body = body_at_depth_one(text, m.end())
        if re.search(r"\bthrow\b", body):
            continue
        patterns = [ASSIGN] + ([DECL] if wide else [])
        for pat in patterns:
            for am in pat.finditer(body):
                name = am.group(1).split("[")[0].split(".")[-1]
                if name not in cond:
                    continue
                line = text[: m.start() + am.start()].count("\n") + 1
                hits.append((line, " ".join(cond.split())[:72], f"{am.group(1)} = {am.group(2)}"))
    return sorted(set(hits))


def main(root, wide):
    sites = files = 0
    for p in sorted(pathlib.Path(root).rglob("*.java")):
        if "/test/" in str(p):
            continue
        hits = scan(p, wide)
        if hits:
            files += 1
            sites += len(hits)
            print(str(p))
            for line, cond, asg in hits:
                print(f"  :{line}   if ({cond})   ->   {asg}")
    label = "WIDE" if wide else "STRICT"
    print(f"\n{label}: {sites} sites in {files} files")


if __name__ == "__main__":
    main(sys.argv[1], "--wide" in sys.argv)
