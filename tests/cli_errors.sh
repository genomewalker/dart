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
removed=(
    --em-streaming --em-max-memory -m --em-null-fdr --coverage-em --coverage-em-iters
    --depth-gate-lo --depth-gate-hi --w-min
    --aln-min-identity --aln-min-bits --aln-max-evalue
    --analysis-prefix --analysis-proteins --analysis-blast8 --analysis-categories
    --blast8-unique --combined-output --corrected --protein-filtered
    --max-dist --threshold --use-identity --no-identity
    --em-iters --em-lambda --em-tol --em-min-margin --em-max-neff --em-audit
    --em-keep-llr --em-reject-llr --em-calib-temperature --preset
)
for f in "${removed[@]}"; do
    expect 1 "was removed" damage-annotate "$f" 1
    expect 1 "was removed" damage-annotate "$f"
done
expect 0 "Usage"                  damage-annotate --no-coverage-em --help
expect 0 "Usage"                  damage-annotate --em --help
expect 0 "Usage"                  damage-annotate --help-all

"$dart" damage-annotate --help >help.txt
"$dart" damage-annotate --help-all >help-all.txt
advanced=(
    --sites --refine-damage --auto-prior-ancient --prior-ancient --d-max --lambda
    --emit-multimap --multimap-k --multimap-floor
    --auto-calibrate-spurious --min-breadth --min-depth --min-reads
    --map --functional-summary --anvio-ko --annotation-source
    --fit-null-from --null-out --em-null-tau --em-checkpoint --resume
    --ancient-threshold --modern-threshold --terminal-threshold --auth-alpha
    --modern-damage-rate --eps0 --per-read-damage-threshold --min-damage-sites
    --site-cap --w0 --em-min-prob --min-positional-score --min-terminal-ratio
    --min-distinct-starts --min-informative-reads --max-single-molecule-frac
    --no-em --em-use-damage-mapping
)
for f in "${advanced[@]}"; do
    if grep -Eq -- "^[[:space:]]+$f([[:space:],]|$)" help.txt; then
        echo "FAIL: default help lists $f"; fail=1
    fi
    if ! grep -Eq -- "(^[[:space:]]+|, )$f([[:space:],]|$)" help-all.txt; then
        echo "FAIL: --help-all is missing $f"; fail=1
    fi
done
for f in "${removed[@]}" --em --no-coverage-em; do
    if grep -Eq -- "(^[[:space:]]+|, )$f([[:space:],]|$)" help.txt help-all.txt; then
        echo "FAIL: help lists removed or hidden flag $f"; fail=1
    fi
done
for f in --emi --damage-index --library-type --no-null --em-null-model --em-calib-model \
    -o --protein-summary --gene-summary --dcx --library-label --em-tape-dir --threads -v; do
    if ! grep -Eq -- "^[[:space:]]+$f([[:space:],]|$)" help.txt; then
        echo "FAIL: default help is missing $f"; fail=1
    fi
done
if [ "$(grep -Ec -- '^[[:space:]]+(-|--)' help.txt)" -ne 14 ]; then
    echo "FAIL: default help must list 14 options"; fail=1
fi

[ "$fail" -eq 0 ] && echo "cli_errors: PASS"
exit "$fail"
