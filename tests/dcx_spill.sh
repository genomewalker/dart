#!/usr/bin/env bash
# damage-annotate --dcx with a forced spill. A 64 KiB DART_DCX_SORT_BYTES makes the external sort
# write runs to disk and merge them; the .dcx, .dcx.auth, .dcx.authstats and the per-read table must
# equal the in-memory sort byte for byte. An injected failure on run 0 proves the spill path ran:
# it must exit 1 with the injected message and leave no run directory behind.
# Usage: dcx_spill.sh <dart binary> <examples/tutorial>. Runs need node-local disk: TMPDIR is used,
# and the test is skipped (77) when TMPDIR is unset, /tmp, or on NFS/tmpfs.
set -u
dart=$1 tut=$2
base=${TMPDIR:-}
case "$base" in ""|/tmp|/tmp/*) echo "SKIP: set TMPDIR to node-local disk (not /tmp)"; exit 77 ;; esac
work=$(mktemp -d "$base/dcx_spill.XXXXXX")
trap 'rm -rf "$work"' EXIT
cd "$work" || exit 1
fail=0
gunzip -c "$tut/expected/hits.emi.gz" > hits.emi
cp "$tut/expected/hits.emi.null" hits.emi.null
mkdir mem spill inject
annotate() {  # tag [env...]
    local tag=$1; shift
    env "$@" "$dart" damage-annotate --emi hits.emi --damage-index "$tut/expected/predictions.agd" \
        --em-tape-dir "$work/$tag" --library-label tutorial -t 4 \
        -o "$tag.tsv" --dcx "$tag.dcx" > "$tag.log" 2>&1
}
annotate mem DART_DCX_TEST_FAIL_RUN=0 || { echo "FAIL: in-memory run exit $?"; tail -3 mem.log; exit 1; }
annotate spill DART_DCX_SORT_BYTES=65536
rc=$?
if [ $rc -ne 0 ] && grep -q 'node-local disk' spill.log; then echo "SKIP: $base is not node-local disk"; exit 77; fi
[ $rc -eq 0 ] || { echo "FAIL: spill run exit $rc"; tail -3 spill.log; exit 1; }
for f in tsv dcx dcx.auth dcx.authstats; do
    if cmp -s "mem.$f" "spill.$f"; then echo "ok: $f identical"; else echo "FAIL: $f differs between in-memory and spilled sort"; fail=1; fi
done
annotate inject DART_DCX_SORT_BYTES=65536 DART_DCX_TEST_FAIL_RUN=0
rc=$?
if [ $rc -ne 1 ] || ! grep -q 'injected run failure' inject.log; then
    echo "FAIL: injected run failure gave exit $rc"; tail -3 inject.log; fail=1
else
    echo "ok: injected failure on run 0 exits 1 (the small budget spilled)"
fi
if ls -d "$work"/inject/dart_dcx_* >/dev/null 2>&1 || ls -d "$work"/spill/dart_dcx_* >/dev/null 2>&1; then
    echo "FAIL: a dart_dcx_* run directory was left behind"; fail=1
else
    echo "ok: run directories removed"
fi
[ $fail -eq 0 ] && echo PASSED
exit $fail
