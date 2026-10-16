# `dart ref-stats` — per-reference coverage stats subcommand (spec)

**Status as of 2026-10-16.** Implemented: `--hits` and `--emi` (EMI v7) inputs, `--threads`, untrimmed
breadth (default) and `--trim`. Validated on the r2s panel and end to end (see "Validation results" below).
Not yet implemented: `--dcx/--ref-lengths`, `depth_std`/`depth_evenness`.
**Alignment length = target span `tend-tstart+1`** in both inputs. It drives `depth_mean` (target residues /
covered) and the `--trim` length. Two things rule out mmseqs `alnlen` here: EMI `aln_len` is not it
(hits2emi-stream stores max(q span, t span)+1, `mmseqs_db.hpp:76`), and damage-annotate reads EMI `aln_len`
as on-target length, so the writer stays unchanged. The lines below that mention `alnlen`/`Σalnlen` are the
original 2026-08-10 spec and are superseded by this paragraph.

**Spec written 2026-08-10.** Motivation: `convergence.env_coverage` breadth is wrong because the convergence
side recomputed coverage from the tlen-less `.dcx` (`condnull/cov/cov_lib.py` → breadth = span **count**).
DART already computes the correct per-reference stats (`src/reference_stats.cpp`, `ReferenceStats`:
`breadth = nonzero_positions / ref_length`), but only inside `annotate` behind `--protein-summary` etc.
Make it a **first-class subcommand** so `env_coverage` is loaded from DART's own numbers — no downstream recompute.

## PER-SAMPLE (one library at a time — NOT pooled)
Coverage stats are computed **per sample/library independently**. Each library is a separate ancient
observation; unioning breadth across the 22 libraries would conflate independent samples. Run the subcommand
once per library → one stats TSV per sample, tagged with `sample_id`. This also makes it fast and
embarrassingly parallel (bounded memory per run; the existing `dist.sh` per-lib-on-chaos-node model), with
**no cross-lib union**.

## CLI (single sample per invocation)
```
dart ref-stats --hits <one_lib.tsv> --sample <LIB> --out <LIB>.cov.tsv [--threads N] [--min-breadth 0]
# or, compact:
dart ref-stats --dcx <one_lib.dcx> --sample <LIB> --ref-lengths <acc\tlen.tsv> --out <LIB>.cov.tsv [--threads N]
```
Driver runs it 22× (one per library, round-robin across chaos nodes), each writing its own `<LIB>.cov.tsv`.
- `--hits`: mmseqs/DART search TSV in the LOCKED 16-col FMT
  `query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln`
  → only cols used: target(1), fident(2), alnlen(3), tstart(8), tend(9), **tlen(13)**, qlen(12). tlen is inline.
- `--dcx + --ref-lengths`: compact path (dcx ~75 MB/lib vs hits ~3.5 GB/lib) — but the `.dcx` **drops tlen**
  (`ga` = query_id,target_id,gamma,…,tstart,tend,…,fident), so lengths must come from `--ref-lengths`
  (build once from `reference_protein` in the ref DB). Prefer this for speed at genome scale.

## Output (reuse `write_reference_stats_tsv`, add `sample_id` + `n_alns`)
`sample_id  ref_id  ref_length  n_reads  n_alns  breadth  depth_mean  depth_std  depth_evenness  avg_identity  avg_alignment_length`
(`breadth` ∈ [0,1] **per sample**; `avg_identity` = mean fident; `depth_mean` = covered-weighted).

## Fast/efficient implementation (do NOT reuse the annotate accumulator as-is)
The existing `CoverageAccumulator` keeps a `std::vector<float>` per position **and** the collector serializes
on **one global mutex** — fine for annotate's per-target subset, poor for a whole library. The subcommand
needs a lean, sharded accumulator, run **per library**:
- **Streaming, mmap + branchless field split**; parse only the needed columns, no per-row allocation.
- **Per-target state without a float depth vector:** a `tlen`-bit **bitset** for covered positions (union of
  `[tstart,tend)`), plus running `sum_fident, sum_alnlen, n_reads, n_alns, min/max start`. breadth =
  popcount/tlen. depth_mean = Σalnlen / covered; skip depth_std/evenness unless a uint16 depth vector is
  requested. **Only this library's targets are resident → bounded memory** (a single lib touches far fewer
  than the genome-wide 47.76 M refs).
- **Target-sharded, lock-free:** hash `target_id` → S shards; one worker per shard owns its map (no mutex).
- **Finalize** per shard → merge → write `<LIB>.cov.tsv`. Optional `--min-breadth` filter.
- **No cross-lib union** — each library is finalized on its own; pooling (if ever wanted) is a separate
  downstream step, never the default.

## Wiring
- `src/cli/cmd_ref_stats.cpp` (new); register via `dart::cli::SubcommandRegistry::instance().register_command(...)`
  (mirror `cmd_predict` in `main.cpp`); add to `CMakeLists.txt` source list (~line 233–243).
- Reuse `ReferenceStats` + `write_reference_stats_tsv` (extend header with `n_alns`).

## Then: load env_coverage PER SAMPLE
`merge_cov.py`/`cov_lib.py` (the pooled `.dcx` recompute) are retired. Run `dart ref-stats` per library →
22 `<LIB>.cov.tsv`, then load **per-sample** rows:
`convergence.env_coverage(release_id, sample_id, env_protein, ref_length, breadth, depth, avg_id, n_alns)`
— note the added `sample_id` and `ref_length`; breadth is already 0–1 **per sample**. The website shows
per-sample coverage for each supporting protein (and can aggregate across samples explicitly if desired, e.g.
max/mean breadth — an explicit reducer, never an implicit union). This is a **schema change** from the current
pooled `env_coverage` (no `sample_id`).

## Validation
- breadth ∈ [0,1] for all rows; spot-check vs `public.r2s_coverage` (400-panel) which already has correct breadth.
- F2D9T0 etc. sanity: breadth = covered/ref_length, not a residue count.

## Validation results (2026-10-16)
Input: the alignments behind `public.r2s_coverage` sample `kapk_r2s_transporters_v1`. That sample is read2Struct
`map-transporters-new/112509/112509_mapping.db`, table `mm_data_filt` (513,455 alns, 904 targets), exported as the
16-col FMT ordered by query.
- Untrimmed `ref-stats` equals an independent SQL union of `[tstart, tend]` on all 904 targets. `n_reads` equals
  the distinct query count on all 904.
- Against `r2s_coverage` (400 proteins), `n_alns` is identical on 400/400 and `avg_identity` is within 1e-4.
  Breadth, compared as covered-position count = round(breadth × ref_length):

  | mode | exact vs r2s_coverage | rest |
  |---|---|---|
  | untrimmed (default) | 0/400 | higher on all; mean +0.12 breadth, Pearson 0.94 |
  | `--trim` | 50/400 | +1 to +14 positions, never lower |

- `r2s_coverage.breadth` is x-filter's `*_cov-stats.tsv.gz`, loaded unchanged by read2Struct (`map.py:506`). It
  differs from ref-stats in two ways:
  1. **Trim** (`x_filter/core_processing.py:652`, on by default, off with `--no-trim`). It zeroes
     `ceil(mean alnLength) + 10` positions at each end but still divides by full length, and is skipped when the
     trim reaches the length. `ref-stats --trim` mirrors this.
  2. **Off-by-one interval** (`x_filter/core_processing.py:498-501`). It uses the 1-based `subjectStart` as a
     0-based index and makes the end exclusive, clamped to `L-1`. This drops the first residue of every alignment
     and the protein's last residue. `ref-stats` does **not** mimic this.
- The trim plus x-filter's interval reproduces `r2s_coverage.breadth` exactly on 400/400 (SQL). So x-filter-era
  breadth runs slightly low, and `ref-stats --trim` is the corrected equivalent.
- Earlier check on the KapK `gold-1000000.tsv` subset (23.1M rows → 10.0M targets, 35 s, 2.1 GB RSS) showed
  breadth/identity in [0,1]. Its breadth is not comparable to the panel (1–4 alns/protein vs 158–1862).
- Outputs are in `/projects/caeg/scratch/kbd606/tmp/dart_refstats_20261016/`.
- The panel numbers above were taken before the switch to target span (mean `alnLength` was used as the
  `--trim` length). Exact `--trim` matches with target span are still 50/400.

## End-to-end validation (2026-10-16)
Pipeline: reads DB `q100000` → existing mmseqs result `res100000` vs `112509-db`
(`faro-ancient-pilots/snipr/kapk-q000-subsets-20261003`). From it, `convertalis` gives `gold-100000.tsv`, and
`hits2emi-stream` gives a v7 `.emi`. `dart predict` is upstream of the result DB, so it is not covered.
Driver: `dart_refstats_20261016/validate_e2e.sh`; summary in `e2e/summary.tsv`.
- **Writer:** `emi-verify` PASS. Rows 2,232,317 = 2,232,317, with 0 multiset differences on
  (query, target, tstart, tend, tlen, qaln, taln) in both directions. Identity max |Δ| 2e-5.
- **ref-stats vs an independent DuckDB oracle over the TSV:** `--hits`/`--emi` × plain/`--trim`, 1,957,244
  targets each. 0 mismatches in ref_length, n_reads, n_alns, breadth, depth_mean and avg_identity. Max |Δ|
  equals the half-unit of output precision.
- **`--hits` vs `--emi`:** identical except avg_identity on 867 targets, at ±1 in the 4th decimal (TSV fident
  has 3 decimals; EMI quantizes to 1/65535).
- **Full library `LV7008956626/hits_v7.emi`:** 3.92 B alignments → 47.57 M targets in 9 min 41 s, 7.6 GB RSS
  at 32 threads. breadth/identity stay in [0,1].
- **Performance (2026-10-16), after adding the sequential prefetch thread and parallel ordered write.**
  Cold timings, NFS, 32 threads, dandycomp07fl:

  | library | size | alignments | open / scan / write | wall | RSS |
  |---|---|---|---|---|---|
  | LV7008956626 (before) | 47 GB | 3.92 B | — | 9 min 41 s (175% CPU) | 7.6 GB |
  | LV7008956626 | 47 GB | 3.92 B | 3.0 / 93.6 / 3.0 s | 1 min 48 s | 7.6 GB |
  | LV7008956628 | 189 GB | 15.39 B | 2.8 / 195.4 / 15.9 s | 3 min 39 s | 11.6 GB |

  The 626 output is byte-identical to the pre-change build's (2.53 GB, 47.57 M rows). The 628 output has
  breadth/identity in [0,1], n_reads ≤ n_alns and no duplicate ref_id. The prefetch pulls the whole file
  through the page cache, so run one library per node.
