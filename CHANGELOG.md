# Changelog

## v2.0.1

Re-pins libtaph to `933877c`, which has the same source as the v2.0.0 pin. Results are unchanged.

## v2.0.0

**Decoy null model.** `predict --fasta-aa X --decoy-every N` also writes `X` with `.decoy`
before the extension: the reversed proteins of every N-th input read, ids prefixed `DECOY_`.
With `--fasta-aa`, N defaults to 10; `--decoy-every 0` disables decoys and an explicit N wins.
`--fasta-aa-masked` alone does not enable decoys. The interval selects input reads, so the
decoy record count includes every ORF of each selected read. Search that file with the real
search's options and pass the hits to
`hits2emi --decoy-hits`, which fits the Gumbel-per-length-bin null and writes it as the binary
sidecar `<emi>.null` with its provenance (dart version, inputs with size and mtime, decoy counts,
`--db-label`). `damage-annotate` loads `<emi>.null` when it exists, after `--em-null-model` and
`--em-calib-model`; `--no-null` is unchanged and is now the only way to run EM without a null.
A null model keeps the EM on its streaming tape: the writer keeps one double per read (the best
null log-density over the read's alignments) and nothing per alignment, so a null run needs the
same `--em-tape-dir` as a `--no-null` run and no more memory.
`dart null` pools decoy hits and sidecars into one null for a batch of samples, and
`dart damage-info FILE.null` prints one. `damage-annotate --fit-null-from/--null-out` are gone; the
TSV sidecar they wrote still loads through `--em-null-model` with a deprecation warning. The fit is
on the raw bit scale and is scaled by the EM's bit-to-nat lambda at load, so the model equals the
old fitter's at the same lambda, except that the old reader dropped the first row of a headerless
decoy file and the new one does not.
Reject bit `1` on the uncalibrated path now comes from the decoy null (`--null-fdr`, or `null_p > --null-p-max`); the old
`gamma_null > --em-null-tau` test now applies only with `--em-calib-model`. The per-read table gains
`null_p` and `null_q` after `gamma_null`. Sidecars are written as version 2, which adds the decoys submitted;
version 1 still loads, with the count unknown.

**Per-bin damaged fraction.** A length bin's `pi_damaged`, `pi_lo` and `pi_hi` used to need
`w_length > 0.5`, which amounted to a coin flip because `w_length` reads 0.500 to 0.502 in
every library. A bin now gets an interval when the damaged-read d_max clears the largest per-bin
δ by 1.96 standard errors.

**`pi_authentic`.** New in each bin of `channels.bulk_damage.bins`: true when the bin has an
interval and `tau_discriminator` is DETECTED. The tau cutoffs (35 and 80 bp) are provisional.

**`w_length` is deprecated.** It is still written but carries no information, since the longest
bin's profile curve is flat and swamps the slope's error. Use `tau_discriminator`. The
preservation terms that read it now use 0, so `qc.preservation.f_coh`, `evidence` and `score`
drop by up to 0.074, 0.27 and 0.18 against earlier builds. No call, d_max source, tau state or
interval count changes. Preservation is diagnostic only and its tier is null.

**tau length only when tau is detected.** `tau_discriminator.point` and
`burial_fingerprint.tau_hat` are null unless `tau_discriminator.state` is DETECTED. Without a
detection the curve is flat from 10 to 400 bp, and the reported length moved by more than 100 bp
between fits that agreed on everything else.

**Faster bulk fit.** The bulk damage fit now uses SQUAREM acceleration and needs a half to a
third of the sweeps on most libraries.

**Capped fits give approximate curves.** tau and the short-bin pi interval read per-bin
profile-likelihood curves, and those are approximate when the main fit stops at the cap. If
`channels.bulk_damage.converged` or a bin's `profile_converged` is false, treat tau and the
short-bin interval as approximate.

**gzip input.** All dart commands decode `.gz` with single-threaded ISA-L, as fqdup does. The
parallel rapidgzip decoder is removed. A `.gz` that ends mid-stream is now an error.

**EM-weighted abundance.** Protein and gene summaries credit each retained candidate to its
own reference with that candidate's EM responsibility (`gamma`). Use `eff_reads` for abundance;
`n_reads` counts candidate assignments and can sum to more than the input read count. Weights
are not renormalized after candidate filtering, so discarded candidates lose their mass.
This applies without `--emit-multimap`; that flag changes per-read output, not summary credit.
Library-wide damage statistics still use one best-hit alignment per read. `--no-em` retains
best-hit assignment. The damaged-read columns (`n_damaged`, `mean_p_damaged`, `frac_damaged`,
`frac_damaged_all`) and the integer authentication counts `auth_binom_k` and `auth_binom_n` are
not abundance columns: each read is counted once, on its best hit, so one damaged read cannot
authenticate several candidate proteins. `max_p_damaged` is taken over every credited read.

**Null-flagged abundance.** When a null rule runs, reads carrying null-reject bit `0x1`
are excluded by default from protein and gene abundance, coverage, assignment counts,
damaged-read counts, authentication counts and summary rows. This applies to both an
explicit `--em-null-model FILE` and an auto-discovered sidecar. `--keep-null-flagged`
restores the previous abundance treatment; without a null model it is a harmless no-op.
Reads carrying only margin, effective-reference
or calibrated-abstain bits keep their existing treatment. The per-read TSV, `.dcx`,
modern-baseline pool, library-wide damage statistics and `null reject:` lines do not change.
The `null drop:` stderr line reports the reads and effective mass actually removed from
abundance; the opt-out instead reports `null drop: off (--keep-null-flagged)`.
This trades sensitivity on short hits for fewer chance-hit assignments.

**`.dcx` sort scratch.** `damage-annotate --dcx` sorts in memory up to 8 GiB and spills sorted
runs past that to `--em-tape-dir` (or `DART_EM_TAPE_DIR`). That directory must be on node-local
disk: NFS, tmpfs and `/tmp` are refused.

**`predict` output no longer depends on the thread count.** The periodic-codon model behind
Pass 2 trained on the first 50,000 eligible reads each thread happened to see, so output varied with
the thread count. The training set is now up to 50,000 reads that are 60 to 200 bp after clipping, every
k-th such read across the whole input (k = their count / 50,000), collected in file
order on one thread; genes are written in input order, and output is identical at any thread
count. The sample costs one extra single-threaded decode of the input,
reported on stderr as `EM training set`. Gene calls therefore differ from earlier versions on
every library; rerun `predict` before comparing or pooling libraries across versions.
