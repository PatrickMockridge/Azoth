#!/usr/bin/env bash
# Run a heavy command inside a memory-capped cgroup scope.
#
# This box has 47 GiB RAM, 2 GiB swap and no OOM guard (systemd-oomd inactive,
# earlyoom absent), so an overrun thrashes the desktop instead of failing one
# process. The kernel only kills after thrash.
#
# Two default allocators are the reason this exists:
#   - a JVM takes 1/4 of RAM as its default heap, measured here at 11.8 GiB;
#     the NeqSim oracle spawns one javac plus ~124 JVMs.
#   - `[profile.release]` sets lto = true, codegen-units = 1, so a release
#     build is a fat-LTO single-codegen-unit link: the peak-memory Rust step.
set -euo pipefail

mem="${AZOTH_GATED_MEMORY_MAX:-12G}"

# Cap the JVM heap unless the caller already chose one.
if [ -z "${JAVA_TOOL_OPTIONS:-}" ] && [[ " $* " == *java* ]]; then
  export JAVA_TOOL_OPTIONS=-Xmx2g
fi

exec systemd-run --user --scope -q -p "MemoryMax=$mem" -p MemorySwapMax=0 "$@"
