# DART: Damage-Aware Read Translation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub release](https://img.shields.io/github/v/release/genomewalker/dart)](https://github.com/genomewalker/dart/releases/latest)
[![CI](https://github.com/genomewalker/dart/actions/workflows/ci.yml/badge.svg)](https://github.com/genomewalker/dart/actions/workflows/ci.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/dart-adna)

DART predicts protein-coding ORFs from ancient DNA reads, writes a damage index, and annotates protein-search alignments. Gene prediction, reference assignment and protein authentication are separate steps.

> Status as of 2026-10-22: checked against the v2.0.0 release branch (libtaph `c3b3a5d`); a null model now runs the EM on the tape path. The decoy null now rejects reads by their chance probability (`null_p`, `--null-p-max`, sidecar version 2); the 0.01 default is under review.

The current main branch uses EMI v7. Damage-aware prediction requires an fqdup library damage summary or explicit self-profiling. `--em-use-damage-mapping` controls the joint damage E-step; protein authentication uses library-level damage information and assigned read evidence.

## Workflow

1. Run `predict` to write GFF3, protein FASTA and an AGD damage index, plus reversed decoy
   peptides with `--decoy-every 10`.
2. Search the proteins and the decoys with MMseqs2; `hits2emi` converts the hits to EMI v7 and
   fits the null model from the decoy hits into `<emi>.null`.
3. Run `damage-annotate` for read assignments, damage scores and protein summaries.
4. Optionally run `ref-stats` for coverage of the input library.

<img src="docs/agp_pipeline.svg" width="700" alt="DART prediction, search and authentication workflow">

## Installation

Requires a C++20 compiler (GCC 10+ or Clang 12+), CMake 3.18+, zlib, zstd, ISA-L, NASM and
pkg-config. With conda:

```bash
conda env create -f environment.yml
conda activate dart
```

Or with system packages: `apt install zlib1g-dev libzstd-dev libisal-dev nasm pkg-config`
(macOS: `brew install zstd isa-l nasm pkg-config`). Then:

```bash
git clone https://github.com/genomewalker/dart.git
cd dart
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

DART's damage estimates come from [libtaph](https://github.com/genomewalker/libtaph), pinned to an
exact commit. If a local checkout exists at `$REPOS/libtaph` (or `-DLIBTAPH_LOCAL_PATH=...`), it
must be at the pinned commit and clean, otherwise configure fails. With no local checkout, CMake
fetches the pinned commit from GitHub. `dart --version` prints the libtaph SHA.

---

## Quick start

The bundled tutorial uses synthetic reads and matching reference proteins:

```bash
cd examples/tutorial
DART=../../build/dart DB=./ref_proteins.faa.gz bash run_tutorial.sh
```

For your own inputs you need one thing the tutorial skips: a damage profile measured with
[fqdup](https://github.com/genomewalker/fqdup) on the reads before deduplication. The null model
that keeps reads from proteins missing from your database off the nearest homolog is built into
the pipeline: `predict --decoy-every 10` writes reversed decoy peptides next to `--fasta-aa`, you
search them with the same MMseqs2 options, `hits2emi --decoy-hits` fits the null into
`hits.emi.null`, and `damage-annotate` finds that file by itself.

```bash
mkdir -p out
# 0) Damage profile from the merged reads, before deduplication
fqdup profile -i merged.fq.gz --damage-summary out/library.tdmg

# 1) Predict genes + damage index. --decoy-every 10 also writes out/proteins.decoy.faa:
#    the reversed proteins of every 10th read, ids prefixed DECOY_
dart predict -i reads.fq.gz -o out/predictions.gff --fasta-aa out/proteins.faa \
  --damage-index out/predictions.agd --library-profile out/library.tdmg --decoy-every 10 -t 16

# 2) Search against reference proteins. VTML20 suits short aDNA peptides; the matrix
#    ships in examples/tutorial/data/ of the DART source tree
SEARCH=(--search-type 1 --min-length 12 -e 10.0 --min-seq-id 0.86 -c 0.65 --cov-mode 2
  --sub-mat examples/tutorial/data/VTML20.out --seed-sub-mat examples/tutorial/data/VTML20.out
  -s 2 -k 6 --spaced-kmer-pattern 11011101
  --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln"
  --threads 16)
mmseqs easy-search out/proteins.faa ref_proteins.faa.gz out/hits.tsv out/tmp "${SEARCH[@]}"

# 3) Search the decoy peptides with the same options (null model input)
mmseqs easy-search out/proteins.decoy.faa ref_proteins.faa.gz out/decoy_hits.tsv out/tmp "${SEARCH[@]}"

# 4) Build the EMI v7 index and fit the null from the decoy hits into out/hits.emi.null
dart hits2emi -i out/hits.tsv -o out/hits.emi --damage-index out/predictions.agd \
  --decoy-hits out/decoy_hits.tsv --decoy-fasta out/proteins.decoy.faa -t 16

# 5) EM assignment + damage annotation (picks up out/hits.emi.null by itself)
dart damage-annotate --emi out/hits.emi --damage-index out/predictions.agd \
  -o out/annotated.tsv --gene-summary out/gene_summary.tsv \
  --auto-calibrate-spurious -t 16
```

One null per database and search configuration is enough for a study: run `--decoy-every 10`
and the decoy search on a few representative samples, pool them with
`dart null -o study.null --db-label mydb sampleA/decoy_hits.tsv sampleB/decoy_hits.tsv ...`,
and pass `--em-null-model study.null` to every `damage-annotate`; the other samples need no
decoys.
`dart damage-info hits.emi.null` prints the fitted bins and the provenance (dart version, inputs,
decoy counts). Refit when the database, the search options, the library type or the read length
range changes.

The null rejects reads by chance probability. For each read's best hit, `null_p` is the
probability that a decoy query gets its best hit in the same aligned-length bin with at least
that bit score: `w_b · (1 − exp(−exp(−(bits − mu_b)/beta_b)))`, where `w_b` is the share of submitted
decoys whose best hit landed in bin `b`. A read with `null_p > --null-p-max` (default 0.01) gets
reject bit `1` and is ambiguous: it leaves the authentication evidence but stays in the
modern-baseline pool, so `q_modern` and the other reads' posteriors do not move. Give the
decoy count with `--decoy-fasta` (or `--decoys-submitted N`; `dart null` takes the latter); without
it `w_b` is per decoy with a hit, which overstates chance. At the end damage-annotate prints one
`null reject:` line per bin with the reads, the share flagged and `S*`, the bit score below which
that bin's reads are flagged (`NA` where `w_b ≤ --null-p-max`, so the bin never flags). Check these
lines: on short peptides the 7-12 aa bins can flag most reads.

The workflow and output interpretation are on the wiki:
[Quick Start Tutorial](https://github.com/genomewalker/dart/wiki/Quick-Start-Tutorial).
Use `--fasta-aa-masked` instead of `--fasta-aa` to replace damage-induced stop codons with `X` so
MMseqs2 can align across them.

### Required choices

| Command | You must pass | Why |
|---------|---------------|-----|
| `predict` | `--library-profile FILE` (binary TDM1 sidecar) or `--allow-self-profile` when damage detection is on | Self-profiling is unreliable on deduplicated libraries |
| `damage-annotate` | a null model when EM is on: `<emi>.null` next to the EMI (written by `hits2emi --decoy-hits`), `--em-null-model FILE`, `--em-calib-model FILE`, or `--no-null` | Without a null model, reads from proteins absent from the database are force-mapped to the nearest homolog |
| `damage-annotate` (EM on) | `--em-tape-dir DIR` or `DART_EM_TAPE_DIR` | Without `--em-calib-model` and `--em-use-damage-mapping`, EM input uses an `O_DIRECT` tape, with or without a null model (the null adds one double per read); supply node-local disk supporting it. Only the calibrated and damage-in-mapping configurations run EM in memory. |

The tutorial opts into self-profiling for its bundled inputs. Its reference database is the 13,681
proteins the reads come from, too few for chance hits: the decoy search yields a handful of hits,
no null bin reaches two decoy queries, and `hits2emi` and `damage-annotate` both say so and run EM
with the null coordinate off. A real database fills the bins. `predict --no-damage` does not require
a profile; `damage-annotate --no-em` bypasses the EM null-model choice.

### Damage profile only

```bash
dart profile -i merged.fq.gz -o profile.json --damage-summary library.tdmg -t 16
```

`dart profile` and `fqdup profile` run the same profiler in libtaph and write byte-identical
JSON and `.tdmg` for the same reads and options. The `.tdmg` is what `predict --library-profile`
reads.

### Changes in v2.0.0

**Decoy null model.** `predict --fasta-aa X --decoy-every N` also writes `X` with `.decoy`
before the extension: the reversed proteins of every N-th read, ids prefixed `DECOY_` (off by
default, since not every predict run feeds a search; 10 is a good value). Search that file with the real search's options and pass the hits to
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
Reject bit `1` on the uncalibrated path now comes from `null_p > --null-p-max`; the old
`gamma_null > --em-null-tau` test could not fire there (the unnormalised real term outscores the
Gumbel density by about 29 nats) and stays only for `--em-calib-model`. The per-read table gains
`null_p` after `gamma_null`. Sidecars are written as version 2, which adds the decoys submitted;
version 1 still loads, with the count unknown.

**Per-bin damaged fraction.** A length bin's `pi_damaged`, `pi_lo` and `pi_hi` used to need
`w_length > 0.5`, which amounted to a coin flip because `w_length` reads 0.500 to 0.502 in
every library. A bin now gets an interval when the damaged-read d_max clears the largest per-bin
δ by 1.96 standard errors. In our test set one library gained intervals and its library-level
`pi` narrowed from [0.35, 0.80] to [0.55, 0.59]. No library lost an interval it had before.

**`pi_authentic`.** New in each bin of `channels.bulk_damage.bins`: true when the bin has an
interval and `tau_discriminator` is DETECTED. The tau cutoffs (35 and 80 bp) were set on three
labelled libraries, so the flag is provisional. Check tau on libraries whose damage status you
already know before relying on it.

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
third of the sweeps on most libraries. Of the 7 test runs that used to stop at the 500-sweep cap,
6 now converge, and no final likelihood went down. On synthetic libraries built to mimic a
pervasive artifact, SQUAREM still stopped at the cap in 8 of 12, so expect some capped fits on
libraries like that.

**Capped fits give approximate curves.** tau and the short-bin pi interval read per-bin
profile-likelihood curves, and those are approximate when the main fit stops at the cap. One
test library still does; its likelihood is within 1e-5 of a 20,000-sweep reference. If
`channels.bulk_damage.converged` or a bin's `profile_converged` is false, treat tau and the
short-bin interval as approximate.

**gzip input.** All dart commands decode `.gz` with single-threaded ISA-L, as fqdup does. The
parallel rapidgzip decoder is gone: fqdup dropped it after it silently truncated a large input
on network storage. A `.gz` that ends mid-stream is now an error. `dart profile` peak memory on
a 2M-read library fell from 1.21 GB to 0.53 GB at `-t 32`.

**`.dcx` sort scratch.** `damage-annotate --dcx` sorts in memory up to 8 GiB and spills sorted
runs past that to `--em-tape-dir` (or `DART_EM_TAPE_DIR`). That directory must be on node-local
disk: NFS, tmpfs and `/tmp` are refused. The `dcx_streaming` test spills under `$TMPDIR`, so run
`ctest` with `TMPDIR` on node-local disk.

**`predict` output no longer depends on the thread count.** The periodic-codon model behind
Pass 2 trained on the first 50,000 eligible reads each thread happened to see, so two runs at
`-t 32` differed in about 0.05% of gene calls and in output order, and `-t 1` differed from
both. The training set is now up to 50,000 reads that are 60 to 200 bp after clipping, every
k-th such read across the whole input (k = their count / 50,000), collected in file
order on one thread; genes are written in input order. Repeated runs at `-t 32` and a run at
`-t 1` are byte-identical. The sample costs one extra single-threaded decode of the input,
reported on stderr as `EM training set`. Gene calls therefore differ from earlier versions on
every library; rerun `predict` before comparing or pooling libraries across versions.

---

## Commands

| Command | Purpose |
|---------|---------|
| `predict` | Gene prediction + damage index (`.agd`) from reads |
| `profile` (`damage`) | Reference-free damage profile (JSON) |
| `hits2emi` | MMseqs2 16-column TSV → EMI v7; `--decoy-hits` fits the null model into `<emi>.null` |
| `null` | Fit or pool the decoy null model (`.null`) shared by the runs of one database |
| `hits2emi-stream` | MMseqs2 result DBs → one EMI v7, streaming, for production-size searches |
| `emi-merge` | Union EMI v7 files over the same reads for joint EM across databases |
| `emi-verify` | Check EMI v7 invariants before annotation |
| `damage-annotate` | EM read assignment, per-read damage posteriors, per-protein authentication |
| `ref-stats` | Per-reference breadth, depth and identity for one sample (`--trim` for x-filter-comparable breadth) |
| `damage-info` | Print the damage profile stored in an `.agd`, or a `.null` null model with its provenance |
| `simulate` | Synthetic aDNA reads (Briggs model) with ground truth |
| `eval` | Score `predict` against `simulate` truth |
| `fit-calib` | Fit the coherent-LR EM calibration model from a labelled search |
| `damage-profile` | Terminal nucleotide profiles of protein-mapped reads |
| `validate` | Validate predictions against aMGSIM ground truth |

Command flags and known help-text discrepancies are documented in the
[Command Reference](https://github.com/genomewalker/dart/wiki/Command-Reference).

### Per-sample coverage

```bash
dart ref-stats --emi SAMPLE.emi --sample SAMPLE --out SAMPLE.cov.tsv --trim --threads 32
```

One row per reference: `sample_id ref_id ref_length n_reads n_alns breadth depth_mean avg_identity`.
`breadth` is covered positions / reference length (0–1). `depth_mean` is aligned target residues / covered positions. Run `ref-stats` separately for each library; measure memory and I/O requirements on your own inputs.

### Production-scale EMI

```bash
dart hits2emi-stream targetDB db1.emi 32 manifest.txt     # manifest: "queryDB resultDB" per line
dart emi-merge joint.emi 32 db1.emi db2.emi
dart emi-verify joint.emi                                  # must print VERIFY: PASS
```

---

## Interpreting results

- **Per read** (`-o`): `posterior` combines terminal (`p_read`), substitution and identity evidence;
  `tier` says which evidence is present. EM columns (`gamma`, `em_keep`, `em_margin`,
  `em_ambiguous`, `em_reject_reason`) say how confidently the read was assigned; ambiguous reads
  do not count as damage evidence.
- **Per protein** (`--protein-summary`): the ancient call is `auth_binom_call`, a binomial test of
  damaged vs modern-rate reads. Trust it when `protein_auth_tier = detected`. A
  `library_not_authenticatable` warning means the library has no usable damage rate, so no protein
  can be called ancient.

Full guide: [Interpreting Results](https://github.com/genomewalker/dart/wiki/Interpreting-Results).
Column schemas: [Output Formats](https://github.com/genomewalker/dart/wiki/Output-Formats).

---

## Benchmarks

The [benchmark directory](benchmark/) contains evaluation scripts and retained results. The [wiki benchmark page](https://github.com/genomewalker/dart/wiki/Benchmarks) describes how to report a versioned comparison. Earlier release measurements are not performance or accuracy guarantees for the current model.

## Citation

If you use DART in your research, please cite:

> **Two-million-year-old microbial communities from the Kap København Formation in North Greenland**
> Fernandez-Guerra A, Wörmer L, Borrel G, et al.
> *bioRxiv* (2025)
> DOI: [10.1101/2023.06.10.544454](https://doi.org/10.1101/2023.06.10.544454)

```bibtex
@article{fernandezguerra2025,
  title={Two-million-year-old microbial communities from the Kap København Formation in North Greenland},
  author={Fernandez-Guerra, Antonio and W{\"o}rmer, Lars and Borrel, Guillaume and Delmont, Tom O and Elberling, Bo and Elvert, Marcus and Eren, A Murat and Gribaldo, Simonetta and Henriksen, Rasmus Amund and Hinrichs, Kai-Uwe and Jochheim, Annika and Korneliussen, Thorfinn S and Krupovic, Mart and Larsen, Nicolaj K and Perez-Laso, Rafael and Pedersen, Mikkel Winther and Pedersen, Vivi K and Ruter, Anthony H and Sand, Karina K and Sikora, Martin and Steinegger, Martin and Veseli, Iva and Wang, Yucheng and Zhao, Lei and {\v{Z}}ure, Marina and Kj{\ae}r, Kurt H and Willerslev, Eske},
  journal={bioRxiv},
  year={2025},
  doi={10.1101/2023.06.10.544454}
}
```

## References

- **Briggs et al. 2007.** Patterns of damage in genomic DNA sequences from a Neandertal. *PNAS* 104(37):14616–14621. [doi:10.1073/pnas.0704665104](https://doi.org/10.1073/pnas.0704665104)
- **Dempster et al. 1977.** Maximum likelihood from incomplete data via the EM algorithm. *J Royal Statistical Society B* 39(1):1–22.
- **Jónsson et al. 2013.** mapDamage2.0. *Bioinformatics* 29(13):1682–1684. [doi:10.1093/bioinformatics/btt193](https://doi.org/10.1093/bioinformatics/btt193)
- **Lindahl 1993.** Instability and decay of the primary structure of DNA. *Nature* 362:709–715. [doi:10.1038/362709a0](https://doi.org/10.1038/362709a0)
- **Michelsen et al. 2022.** metaDMG. *bioRxiv*. [doi:10.1101/2022.12.06.519264](https://doi.org/10.1101/2022.12.06.519264)
- **Müller & Vingron 2000.** Modeling amino acid replacement. *J Computational Biology* 7(6):761–776.
- **Steinegger & Söding 2017.** MMseqs2. *Nature Biotechnology* 35(11):1026–1028. [doi:10.1038/nbt.3988](https://doi.org/10.1038/nbt.3988)
- **Van der Jeugt et al. 2022.** FragGeneScanRs. *BMC Bioinformatics* 23(1):198. [doi:10.1186/s12859-022-04736-5](https://doi.org/10.1186/s12859-022-04736-5)

---

## License

MIT License; see LICENSE file for details.
