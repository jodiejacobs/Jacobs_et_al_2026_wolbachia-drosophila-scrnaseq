# Jacobs et al. 2026 — *Wolbachia*–*Drosophila* scRNA-seq

Snakemake pipeline and analysis scripts for single-cell RNA-seq profiling of *w*Mel *Wolbachia*-infected *Drosophila melanogaster* JW18 cell lines. The pipeline aligns reads against a combined *D. melanogaster* + *w*Mel + 16S reference with kallisto/bustools, filters and QCs the resulting count matrices, integrates samples across sequencing platforms (10x Genomics and PIP-seq) and infection conditions, and runs downstream analyses including cell cycle assignment, NMF gene program discovery, cluster marker/pathway analysis, pseudotime gene importance (tradeSeq), and titer–expression correlation.

## Repository layout

- `cell_culture_system/` — the Snakemake pipeline for the JW18 cell culture experiment.
  - `Snakefile` — main pipeline (rules listed below).
  - `config/config.yaml` — paths, conda env locations, and per-rule SLURM resource settings.
  - `config/samples.csv` — sample sheet (condition, platform, replicate, R1 path, R2 path).
  - `config/env.yml`, `config/envs/scanpy_env.yml`, `config/envs/cyclum_env.yml` — conda/mamba environment specs.
  - `snakevision/` — rule graph visualization (`pipeline_rulegraph.svg`) and the command used to generate it.
- `snakemake_scripts/` — scripts called by the Snakefile, grouped by pipeline stage:
  - `alignment/` — kb-python/kallisto alignment notes and BUS-file postprocessing.
  - `filtering/` — QC filtering of kallisto/bustools AnnData objects.
  - `quality_control/` — QC plots and metrics on aligned/filtered data.
  - `method_comparison/` — integration across 10x/PIP-seq platforms, cell cycle scoring, cluster marker and pathway analysis, pseudotime (tradeSeq), NMF program discovery/annotation, SCEPTIC, and platform concordance validation.
  - `analysis/` — cell cycle (Cyclum) analysis, NMF program annotation, GSEA, Wolbachia titer analysis.
  - `plotting/` — summary statistics and QC histogram plotting.
- `parsing_scripts/` — standalone scripts for extracting/summarizing 16S and *Wolbachia* titer reads from raw data, independent of the main Snakemake run (`analysis/`, `data_processing/`, `extract_by_barcode_list.py`, `notes_for_16S.txt`).
- `blast_db/` — prebuilt BLAST database (16S ribosomal RNA + *Wolbachia*) used to classify captured 16S reads.

## Requirements

- Linux with [mamba](https://github.com/mamba-org/mamba)/conda
- Snakemake ≥ 9.0
- SLURM (the Snakefile is written to submit jobs via `snakemake --executor slurm`; adjust `config.yaml` resource/partition settings for other schedulers)
- kallisto, bustools, kb-python (see `snakemake_scripts/alignment/kb_notes.txt`)
- Two analysis conda environments, referenced in `config.yaml`:
  - `scanpy_env` — scanpy/anndata-based filtering, QC, integration, clustering, NMF, GSEA (spec in `config/envs/scanpy_env.yml` / `config/env.yml`)
  - `cyclum_env` — Cyclum-based cell cycle scoring (spec in `config/envs/cyclum_env.yml`)
- R with `tradeSeq` for pseudotime gene importance (`run_tradeseq.R`, `run_tradeseq_hclust.R`)
- NCBI BLAST+ (for 16S read classification against `blast_db/`)

## Setup

```bash
# Alignment environment
mamba create -n kallisto_env -c bioconda kallisto bustools gffread kb-python

# Analysis environments (edit paths in config.yaml to match)
mamba env create -f cell_culture_system/config/envs/scanpy_env.yml
mamba env create -f cell_culture_system/config/envs/cyclum_env.yml
```

Edit `cell_culture_system/config/config.yaml` to point `scanpy_env`/`cyclum_env` at your environment paths, and `cell_culture_system/config/samples.csv` at your FASTQ locations. Build the combined reference index per `snakemake_scripts/alignment/kb_notes.txt` and point `kallisto_index` / `transcripts_to_genes` at it.

## Sample sheet format

`config/samples.csv` is a headerless, comma-separated file:

```
condition, platform, replicate, R1_path, R2_path
```

`platform` is `10x` or `pipseq`. Conditions in this study include *Wolbachia*-infected (`JW18wMel`) and antibiotic-cured (`JW18DOX`) lines, sampled as untreated controls and across a *Sindbis virus* infection time course (`-SV-D1`, `-D7`, `-D28`, `-D56`).

## Running the pipeline

```bash
mamba activate snakemake  # requires snakemake >= 9.0
cd cell_culture_system

# Dry run
snakemake --executor slurm \
  --default-resources slurm_partition=medium slurm_time="2:00:00" runtime=120 mem_mb=8000 \
  -j 16 -n

# Full run
snakemake --executor slurm \
  --default-resources slurm_partition=medium slurm_time="2:00:00" runtime=120 mem_mb=8000 \
  -j 16
```

Per-rule thread/memory/time/partition settings live in `config.yaml` and override the defaults above.

### Pipeline stages (Snakefile rules)

1. **Alignment** — `map_pipseq`, `align_gene_reads` (kallisto/bustools pseudoalignment against the *D. melanogaster* + *w*Mel + 16S reference)
2. **16S / titer QC** — `calculate_coverage`, `summarize_blast`, `plot_coverage_by_group`, `plot_blast_by_group` (extraction and BLAST classification of captured 16S rRNA reads)
3. **Combine & integrate** — `combine_files_by_condition_platform`, `integrate`, `integrate_uninfected` (merge replicates/platforms into condition-level AnnData objects)
4. **Cell cycle** — `annotate_cell_cycle`, `cell_cycle_analysis`, `cell_cycle_analysis_uninfected`, `cyclum_analysis_uninfected`, `project_to_cell_cycle`
5. **Clustering & markers** — `cluster_marker_pathway`
6. **Gene programs** — `nmf_programs`, `nmf_continuous_var`, `nmf_categorical_var`, `nmf_annotate_programs`
7. **Trajectory / titer modeling** — `run_sceptic`, `pseudotime_gene_importance`, `run_tradeseq`
8. **Validation** — `validate_platform_concordance` (10x vs. PIP-seq concordance)
9. `rule all` ties the full DAG together.

A rendered rule graph is in `cell_culture_system/snakevision/pipeline_rulegraph.svg` (regenerate per the command in `snakevision/readme`).

## 16S / *Wolbachia* titer extraction (standalone)

`parsing_scripts/notes_for_16S.txt` documents the manual workflow for pulling *Wolbachia* 16S reads out of a BUS file (`bustools capture` → `bustools text` → `filter_fastq.py`) for downstream QC/mapping outside the main Snakemake DAG. Titer quantification scripts live in `parsing_scripts/analysis/wolbachia_titer_analysis.py` and its batch wrappers.

## Citation

If you use this pipeline, please cite Jacobs et al. 2026 (citation details to be added on publication).