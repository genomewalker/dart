#!/usr/bin/env bash
# Bad input must give a message and exit 1, never an abort (134) or a silent exit 0.
# Usage: cli_errors.sh <dart binary>
set -u
dart=$1
work=$(mktemp -d "${TMPDIR:-.}/cli_errors.XXXXXX")
trap 'rm -rf "$work"' EXIT
cd "$work" || exit 1
ulimit -c 0
fail=0

# expect <exit> <stderr pattern> <args...>
expect() {
    local want=$1 pat=$2; shift 2
    "$dart" "$@" >out.txt 2>err.txt
    local got=$?
    if [ "$got" -ne "$want" ] || ! grep -q -- "$pat" err.txt out.txt; then
        echo "FAIL: dart $* -> exit $got (want $want, pattern '$pat')"; sed 's/^/    /' err.txt | head -3
        fail=1
    fi
}

expect 1 "cannot open"            emi-verify missing.emi
expect 1 "cannot open"            emi-merge out.emi 1 missing1.emi missing2.emi
expect 1 "cannot open"            damage-annotate --emi missing.emi --no-null
expect 1 "invalid numeric value"  damage-annotate --emi missing.emi --em-min-prob abc
expect 1 "invalid numeric value"  damage-profile --min-reads abc
expect 1 "must be ss, ds or auto" damage-annotate --emi missing.emi --library-type xx
expect 1 "cannot open"            eval --truth missing.tsv missing.gff
expect 1 "cannot open"            validate missing.fq.gz missing.tsv.gz missing.gff
for f in --em-max-memory --em-null-fdr --coverage-em --coverage-em-iters --depth-gate-lo --w-min; do
    expect 1 "was removed" damage-annotate "$f" 1
done
expect 0 "Usage"                  damage-annotate --no-coverage-em --help

if "$dart" damage-annotate --help 2>&1 | grep -E -- "--em-streaming|--em-max-memory|--coverage-em|--em-null-fdr|--depth-gate|--w-min"; then
    echo "FAIL: damage-annotate --help still lists a removed flag"; fail=1
fi

[ "$fail" -eq 0 ] && echo "cli_errors: PASS"
exit "$fail"
