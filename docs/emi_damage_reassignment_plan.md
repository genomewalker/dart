# EMI Damage-Aware Reassignment Plan (Principled)

## Objective
Build a mathematically correct reassignment workflow where damage evidence can change multi-mapping assignments when it is target-specific.

## Core Principle
- Separate per-read prior from per-alignment evidence.
- `p_r` is read-level prior (`P(read r is ancient)`), from AGD.
- `d_rt` is alignment-level damage evidence (`read r` vs `target t`), computed per alignment.
- `d_max/lambda` in EMI header remain sample-level parameters, not per-alignment discriminators.

## Problem in Current Behavior
- EMI `DAMAGE_SCORE` is currently copied from one AGD value per read to all alignments of that read.
- This creates no within-read variation, so damage cancels in EM normalization and cannot affect reassignment.

## Target EM Model
For each read `r` and candidate target `t`:

`log w_rt = log(theta_t) + lambda_b * bit_rt - 0.5*log(L_t) + log( p_r * exp(llA_rt) + (1 - p_r) * exp(llM_rt) )`

Then:

`gamma_rt = softmax_t(log w_rt)`

M-step:

`theta_t ∝ sum_r gamma_rt + (alpha - 1)`

Notes:
- `llA_rt = log P(damage-data_rt | ancient, t)`
- `llM_rt = log P(damage-data_rt | modern, t)`
- Equivalent storage option: `logBF_rt = llA_rt - llM_rt` plus a stable implementation of the mixture term.

## EMI v5 Schema Changes
- Add `DMG_K` (uint16): damage-consistent terminal hits per alignment.
- Add `DMG_M` (uint16): terminal opportunities per alignment.
- Add `DMG_LL_A` (float): log-likelihood under ancient model.
- Add `DMG_LL_M` (float): log-likelihood under modern model.
- Keep existing `DAMAGE_SCORE` as read-level prior if needed for compatibility in summaries.
- Bump EMI version to 5 and document migration rules.

## Implementation Plan
1. EMI format and IO
- [ ] Add new columns and metadata in EMI header/column enum.
- [ ] Update writer path in `hits2emi` to compute and persist per-alignment damage evidence.
- [ ] Update reader scan APIs and column selection to expose the new columns.
- [ ] Add compatibility path for older EMI versions with explicit warning.

2. Damage evidence computation in `hits2emi`
- [ ] Compute per-alignment terminal opportunities and hits from alignment strings.
- [ ] Compute `llA_rt` and `llM_rt` from `d_max/lambda` and baseline mismatch model.
- [ ] Cache reusable per-read/per-position components to avoid repeated hot-loop work.
- [ ] Ensure this stage remains parallel and vectorization-friendly.

3. EM reassignment update
- [ ] Replace current damage gating with unconditional use of alignment-level evidence when present.
- [ ] Implement stable log-space mixture term using `p_r`, `llA_rt`, `llM_rt`.
- [ ] Remove global `pi` update from EM reassignment path (not needed with fixed `p_r` prior).
- [ ] Keep SQUAREM and objective safeguards.

4. Diagnostics and logs
- [ ] Print evidence quality stats: `%alignments with DMG_M>0`, distribution of `logBF_rt`.
- [ ] Print EM diagnostics with human-readable times and objective deltas.
- [ ] Explicitly report whether damage is informative for reassignment and why.

5. Output contract for downstream analysis
- [ ] Keep single primary annotated TSV.
- [ ] Keep optional machine-friendly outputs with deterministic names for workflow engines.
- [ ] Add columns exposing EM posterior (`gamma`), reassignment flag, and damage evidence summary.

6. Validation and QA
- [ ] Unit tests: numerical stability for log-space computations and edge cases (`p_r` near 0/1).
- [ ] Property tests: when `llA_rt == llM_rt` for all `t`, damage must have no reassignment effect.
- [ ] Integration tests: synthetic multi-mapper where one target is damage-consistent, one is not.
- [ ] Regression benchmark on KEGG-scale dataset.

7. Performance goals
- [ ] No repeated string parsing in EM hot loop.
- [ ] Numeric-only EM inputs after EMI build.
- [ ] Maintain parallel scans and per-thread accumulation.
- [ ] Track I/O vs CPU time separately in logs (important for NFS runs).

## Acceptance Criteria
- Multi-mapping reassignment changes when and only when alignment-level damage evidence differs across candidates.
- EM remains numerically stable and converges monotonically with safeguards.
- Runtime is not worse than current pipeline on equivalent hardware, excluding unavoidable storage bottlenecks.
- Outputs remain directly consumable by existing downstream scripts without manual reshaping.

