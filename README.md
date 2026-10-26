# DART: Damage-Aware Read Translation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub release](https://img.shields.io/github/v/release/genomewalker/dart)](https://github.com/genomewalker/dart/releases/latest)
[![CI](https://github.com/genomewalker/dart/actions/workflows/ci.yml/badge.svg)](https://github.com/genomewalker/dart/actions/workflows/ci.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/dart-adna)

DART predicts protein-coding ORFs from ancient DNA reads, writes a damage index, and annotates protein-search alignments. Gene prediction, reference assignment and protein authentication are separate steps.

DART uses the EMI v7 index format. Damage-aware prediction requires an fqdup library damage summary or explicit self-profiling. `--em-use-damage-mapping` controls the joint damage E-step; protein authentication uses library-level damage information and assigned read evidence.

## Workflow

1. Run `predict` to write GFF3, protein FASTA and an AGD damage index, plus reversed decoy
   peptides by default with `--fasta-aa` (`--decoy-every 10`; `0` opts out).
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

Or with system packages: `apt install zlib1g-dev libzstd-dev libisal-dev nasm pkg-config`.
Then:

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

Tests: `ctest --test-dir build`, with `TMPDIR` on a local disk (not NFS or tmpfs).

---

## Quick start

For many samples, [dart-nf](https://github.com/genomewalker/dart-nf) runs the same steps (damage profile, predict, MMseqs2 target and decoy searches, hits2emi, damage-annotate) as a Nextflow pipeline.

The bundled tutorial uses synthetic reads and matching reference proteins:

```bash
cd examples/tutorial
DART=../../build/dart DB=./ref_proteins.faa.gz bash run_tutorial.sh
```

For your own inputs you need one thing the tutorial skips: a damage profile measured with
[fqdup](https://github.com/genomewalker/fqdup) on the reads before deduplication. The null model
that keeps reads from proteins missing from your database off the nearest homolog is built into
the pipeline: `predict --fasta-aa FILE` writes reversed decoy peptides by default (every 10th
input read). Search them with the same MMseqs2 options, `hits2emi --decoy-hits` fits the null into
`hits.emi.null`, and `damage-annotate` finds that file by itself.

```bash
mkdir -p out
# 0) Damage profile from the merged reads, before deduplication
fqdup profile -i merged.fq.gz --damage-summary out/library.tdmg

# 1) Predict genes + damage index. Also writes out/proteins.decoy.faa by default:
#    the reversed proteins of every 10th read, ids prefixed DECOY_
dart predict -i reads.fq.gz -o out/predictions.gff --fasta-aa out/proteins.faa \
  --damage-index out/predictions.agd --library-profile out/library.tdmg -t 16

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
and pass `--em-null-model study.null` to every `damage-annotate`; the other samples can opt out
with `predict --decoy-every 0`.
`dart damage-info hits.emi.null` prints the fitted bins and the provenance (dart version, inputs,
decoy counts). Refit when the database, the search options, the library type or the read length
range changes.

The null rejects reads by chance probability. For each read's best hit, `null_p` is the
probability that a decoy query gets its best hit in the same aligned-length bin with at least
that bit score: `w_b · (1 − exp(−exp(−(bits − mu_b)/beta_b)))`, where `w_b` is the share of submitted
decoys whose best hit landed in bin `b`. With `Q` queries searched, `Q · null_p` is the number of
chance best hits expected at or above that score in the bin. Within each bin damage-annotate ranks
the reads by bit score, divides the expected chance hits at each score by the reads at or above it,
and takes the running minimum from below: that is the read's q-value, the false discovery rate of
keeping it and every read above it, written per read as `null_q` (`NA` without a hit, in an unfit
bin, or under `--null-p-max`), so a stricter cutoff is a filter on one run's table. A read with q above `--null-fdr` (default 0.05) gets reject bit
`1` and is ambiguous: it leaves the authentication evidence but stays in the modern-baseline pool,
so `q_modern` and the other reads' posteriors do not move. `Q` is the record count of
`--damage-index` (one record per predicted protein, the proteins searched); pass
`--null-queries N` when the searched set differs. Every query counts, including those with a true
hit, so the rule is conservative. Give the decoy count with `--decoy-fasta` (or
`--decoys-submitted N`; `dart null` takes the latter). `hits2emi --decoy-hits` requires one
of these options; omitting both is an error. Older sidecars with an unknown count still use
decoys with a hit as the denominator, which overstates chance. At the end damage-annotate prints
one `null reject:` line per bin with
the reads, the share flagged, `chance` (the chance best hits expected in the bin, `Q · w_b`) and the
lowest kept bit score (`NA` when nothing in the bin is flagged). A bin with no decoy hits of its own
never flags. `--null-p-max P` instead rejects a read when `null_p > P`, a fixed per-query chance;
only a bin holding more than `P` of the decoys can then flag, so on a large database a threshold
like 0.05 flags nothing. The two options are exclusive. Without a damage index or `--null-queries`
damage-annotate falls back to `null_p > 0.05` and says so.

The workflow and output interpretation are on the wiki:
[Quick Start Tutorial](https://github.com/genomewalker/dart/wiki/Quick-Start-Tutorial).
Use `--fasta-aa-masked` instead of `--fasta-aa` to replace damage-induced stop codons with `X` so
MMseqs2 can align across them.

### Required choices

| Command | You must pass | Why |
|---------|---------------|-----|
| `predict` | `--library-profile FILE` (binary TDM1 sidecar) or `--allow-self-profile` when damage detection is on | Self-profiling is unreliable on deduplicated libraries |
| `hits2emi --decoy-hits FILE` | `--decoy-fasta FILE` or `--decoys-submitted N` (> 0) | The null denominator counts all decoys searched, including those without hits |
| `damage-annotate` | a null model when EM is on: `<emi>.null` next to the EMI (written by `hits2emi --decoy-hits`), `--em-null-model FILE`, `--em-calib-model FILE`, or `--no-null` | Without a null model, reads from proteins absent from the database are force-mapped to the nearest homolog |
| `damage-annotate` (EM on) | `--em-tape-dir DIR` or `DART_EM_TAPE_DIR` | Without `--em-calib-model` and `--em-use-damage-mapping`, EM input uses an `O_DIRECT` tape, with or without a null model (the null adds one double per read); supply node-local disk supporting it. Only the calibrated and damage-in-mapping configurations run EM in memory. |

The tutorial opts into self-profiling for its bundled inputs. Its reference database is the 13,681
proteins the reads come from, too few for chance hits: the decoy search yields a handful of hits,
no null bin reaches two decoy queries, and `hits2emi` and `damage-annotate` both say so and run EM
with the null coordinate off. A real database fills the bins. `predict --no-damage` does not require
a profile; `damage-annotate --no-em` bypasses the EM null-model choice.

`damage-annotate --lambda F` sets the decay constant explicitly. By default it comes from
the damage index (or the EMI header when no damage index is supplied); `0.3` is the fallback
when absent.

### Damage profile only

```bash
dart profile -i merged.fq.gz -o profile.json --damage-summary library.tdmg -t 16
```

`dart profile` and `fqdup profile` run the same profiler in libtaph and write byte-identical
JSON and `.tdmg` for the same reads and options. The `.tdmg` is what `predict --library-profile`
reads.

Release notes: [CHANGELOG.md](CHANGELOG.md).

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
