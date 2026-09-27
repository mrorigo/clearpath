#!/usr/bin/env bash
# The gates from `docs/SPEC.md` section 10, runnable without CI.
#
# Each is a hard failure, not a warning. The `unsafe` gate is a grep rather than a lint because
# `#![deny(unsafe_code)]` cannot express "unsafe is allowed in exactly one file".
set -euo pipefail
cd "$(dirname "$0")/.."

fail() { echo "FAIL: $*" >&2; exit 1; }
pass() { echo "  ok: $*"; }

echo "==> build"
cargo build --all-targets
pass "build"

echo "==> no_std"
cargo build --no-default-features
pass "no_std (the crate is #![no_std]; --no-default-features drops the std dependency)"

# A real cross-compile proves more, and is worth doing when the target is installed.
for target in thumbv7em-none-eabihf aarch64-unknown-none; do
    if rustup target list --installed 2>/dev/null | grep -qx "$target"; then
        cargo build --no-default-features --target "$target"
        pass "no_std cross-compile to $target"
    else
        echo "  skip: $target not installed (rustup target add $target)"
    fi
done

echo "==> clippy"
cargo clippy --all-targets --all-features -- -D warnings
pass "clippy -D warnings"

echo "==> unsafe containment"
# Only the SIMD module may contain `unsafe`, and only under an explicit allow. The grep is for
# `unsafe` in *code* positions: a doc comment that merely names `unsafe` is not a violation, and
# matching one would make the gate cry wolf.
violations="$(grep -rnE --include='*.rs' '(unsafe[[:space:]]*\{|unsafe[[:space:]]+(fn|impl|trait|extern))' src/ \
    | grep -v '^src/spline/simd.rs:' || true)"
if [ -n "$violations" ]; then
    echo "$violations" >&2
    fail "unsafe outside src/spline/simd.rs"
fi
if [ -f src/spline/simd.rs ] && ! grep -q 'allow(unsafe_code)' src/spline/simd.rs; then
    fail "src/spline/simd.rs has no #![allow(unsafe_code)]"
fi
pass "no unsafe outside src/spline/simd.rs"

echo "==> tests"
cargo test --all-features
pass "tests"

echo "==> tests, repeated"
# The property tests reseed every run, and two of the bugs found during development only
# reproduced on some seeds. One green run is not evidence for either.
for _ in 1 2 3 4 5 6; do
    cargo test --all-features --quiet || fail "a repeated test run failed"
done
pass "6 repeated runs"

echo "==> doc tests"
cargo test --doc
pass "doc tests"

echo
echo "all gates passed"
