# DART: Damage-Aware Read Translation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub release](https://img.shields.io/github/v/release/genomewalker/dart)](https://github.com/genomewalker/dart/releases/latest)
[![CI](https://github.com/genomewalker/dart/actions/workflows/ci.yml/badge.svg)](https://github.com/genomewalker/dart/actions/workflows/ci.yml)

Ancient DNA reads are short and carry deamination damage: C→T near the 5′ end and G→A (or C→T)
near the 3′ end. These changes create false stop codons and wrong amino acids, so standard gene
callers and protein searches lose or misplace ancient reads. DART translates damaged reads into
proteins, assigns each read to the reference protein it most likely came from, and tests whether
each protein carries the damage pattern of ancient DNA.

## Workflow

1. `predict` translates the reads, writing GFF3, protein FASTA, a damage index (`.agd`) and
   reversed decoy peptides used to model chance hits.
2. MMseqs2 searches the proteins and the decoys against a reference database.
3. `hits2emi` indexes the hits (EMI format) and fits the chance-hit model from the decoy hits.
4. `damage-annotate` resolves multi-mapping reads with EM, scores damage per read and tests
   each protein for ancient damage.

<img src="docs/agp_pipeline.svg" width="700" alt="DART prediction, search and authentication workflow">

To run many samples, use [dart-nf](https://github.com/genomewalker/dart-nf), a Nextflow pipeline
for the whole workflow.

## Installation

Requires a C++20 compiler (GCC 10+ or Clang 12+), CMake 3.18+, zlib, zstd, ISA-L, NASM and
pkg-config.

```bash
conda env create -f environment.yml && conda activate dart   # or: apt install zlib1g-dev libzstd-dev libisal-dev nasm pkg-config
git clone https://github.com/genomewalker/dart.git
cd dart
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

CMake fetches [libtaph](https://github.com/genomewalker/libtaph), which provides the damage
estimates, at a pinned commit; `dart --version` prints it. To build against a local checkout,
pass `-DLIBTAPH_LOCAL_PATH=...`; it must be clean and at the pinned commit. Run the tests with
`ctest --test-dir build`, with `TMPDIR` on a local disk.

## Quick start

The bundled tutorial runs the full workflow on synthetic reads:

```bash
cd examples/tutorial
DART=../../build/dart DB=./ref_proteins.faa.gz bash run_tutorial.sh
```

On your own data:

```bash
mkdir -p out
# Damage profile from the merged reads, before deduplication
fqdup profile -i merged.fq.gz --damage-summary out/library.tdmg

# Predict genes; also writes out/proteins.decoy.faa (reversed proteins of every 10th read)
dart predict -i reads.fq.gz -o out/predictions.gff --fasta-aa out/proteins.faa \
  --damage-index out/predictions.agd --library-profile out/library.tdmg -t 16

# Search proteins and decoys with the same options (VTML20 suits short ancient peptides)
SEARCH=(--search-type 1 --min-length 12 -e 10.0 --min-seq-id 0.86 -c 0.65 --cov-mode 2
  --sub-mat examples/tutorial/data/VTML20.out --seed-sub-mat examples/tutorial/data/VTML20.out
  -s 2 -k 6 --spaced-kmer-pattern 11011101
  --format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qaln,taln"
  --threads 16)
mmseqs easy-search out/proteins.faa ref_proteins.faa.gz out/hits.tsv out/tmp "${SEARCH[@]}"
mmseqs easy-search out/proteins.decoy.faa ref_proteins.faa.gz out/decoy_hits.tsv out/tmp "${SEARCH[@]}"

# Index the hits and fit the chance-hit model into out/hits.emi.null
dart hits2emi -i out/hits.tsv -o out/hits.emi --damage-index out/predictions.agd \
  --decoy-hits out/decoy_hits.tsv --decoy-fasta out/proteins.decoy.faa -t 16

# Assign reads, score damage, summarise genes and proteins
dart damage-annotate --emi out/hits.emi --damage-index out/predictions.agd \
  -o out/annotated.tsv --gene-summary out/gene_summary.tsv \
  --protein-summary out/protein_summary.tsv --auto-calibrate-spurious \
  --em-tape-dir /local/disk/tmp -t 16
```

Things to know:

- **Damage profile.** `predict` needs a library profile from `fqdup profile` or `dart profile`,
  measured before deduplication, or `--allow-self-profile`.
- **Chance hits.** Reads whose best hit could be a chance match, judged against the decoy
  hits (`null_q` above `--null-fdr`, default 0.05), are flagged and left out of gene and protein abundance. This also removes some
  true short hits: on a test sample of known content, about half of the true 10–12 aa hits.
  Pass `--keep-null-flagged` to keep them. Details:
  [Spurious-Hit Filters](https://github.com/genomewalker/dart/wiki/Spurious-Hit-Filters).
- **One null per study.** The chance-hit model depends on the database and search options, not
  the sample. Fit it on a few samples, pool them with `dart null -o study.null ...`, and pass
  `--em-null-model study.null` to the rest (which can skip decoys with `--decoy-every 0`).
- **EM scratch.** `damage-annotate` streams the EM through `--em-tape-dir` (or
  `DART_EM_TAPE_DIR`), which must be on local disk, not NFS or tmpfs.

## Commands

| Command | Purpose |
|---------|---------|
| `predict` | Gene prediction and damage index (`.agd`) from reads |
| `profile` | Reference-free damage profile (JSON and `.tdmg`) |
| `hits2emi` | MMseqs2 TSV → EMI index; `--decoy-hits` fits the chance-hit model |
| `hits2emi-stream` | MMseqs2 result databases → EMI, streaming, for large searches |
| `emi-merge` | Merge EMI files over the same reads for joint EM across databases |
| `emi-verify` | Check an EMI before annotation |
| `null` | Fit or pool the chance-hit model for a set of samples |
| `damage-annotate` | EM read assignment, per-read damage and per-protein authentication |
| `ref-stats` | Per-reference breadth, depth and identity for one sample |
| `damage-info` | Print the damage profile in an `.agd` or a `.null` model |
| `simulate` | Synthetic ancient reads (Briggs model) with ground truth |
| `eval` | Score `predict` against `simulate` truth |
| `fit-calib` | Fit the EM calibration model from a labelled search |
| `damage-profile` | Terminal nucleotide profiles of protein-mapped reads |
| `validate` | Validate predictions against aMGSIM ground truth |

Flags and defaults: [Command Reference](https://github.com/genomewalker/dart/wiki/Command-Reference).

## Outputs

- **Per read** (`-o`): `posterior` is the probability that the read is ancient, combining terminal
  damage, damage-consistent substitutions and identity. The EM columns (`gamma`, `em_keep`,
  `em_margin`, `em_reject_reason`) say how confidently the read was assigned.
- **Per gene and protein** (`--gene-summary`, `--protein-summary`): use `eff_reads` for
  abundance. The ancient call is `auth_binom_call`; trust it when `protein_auth_tier` is
  `detected`. A `library_not_authenticatable` warning means the library has no usable damage
  rate, so no protein can be called ancient.

Column definitions: [Output Formats](https://github.com/genomewalker/dart/wiki/Output-Formats)
and [Interpreting Results](https://github.com/genomewalker/dart/wiki/Interpreting-Results).
Model details: [Methods and Model](https://github.com/genomewalker/dart/wiki/Methods-and-Model).
Release notes: [CHANGELOG.md](CHANGELOG.md).

## Benchmarks

Evaluation scripts and results are in [benchmark/](benchmark/); see the
[Benchmarks](https://github.com/genomewalker/dart/wiki/Benchmarks) wiki page.

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

## License

MIT. See [LICENSE](LICENSE).
