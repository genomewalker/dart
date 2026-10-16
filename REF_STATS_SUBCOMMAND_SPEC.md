# `dart ref-stats` — per-reference coverage stats subcommand (spec)

**Status as of 2026-08-10.** Motivation: `convergence.env_coverage` breadth is wrong because the convergence
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
