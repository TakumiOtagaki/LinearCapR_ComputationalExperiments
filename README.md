# LinCapR Experiments (Reproducibility Guide)

This repository packages the scripts, configs, and helper binaries we used to generate the figures in the manuscript comparing **CapR** and **LinCapR**. The goal is to let others rebuild the publication figures from the provided inputs.

## Prerequisites
- Git submodules: `git submodule update --init --recursive` (brings in `CapR/` and `LinearCapR/`).
- Build CapR / LinCapR:
  - CapR: `cd CapR && make` (adjust compiler if needed).
  - LinCapR: `cd LinearCapR && make`.
- Python Environment: `uv sync`

## Data layout
We assume the processed datasets already exist in `data/processed/`:
- bpRNA multi-FASTA: `data/processed/bprna/bpRNA_multifasta.fasta`
- RNAcentral sample: `data/processed/rnacentral/rnacentral_active_20.fasta`
- SARS-CoV-2: `data/processed/sarscov2/NC_045512.RNA.fa`

Outputs from earlier runs (ROC, time/memory, etc.) are referenced under `result/` and `graphs/`. If you have your own runs, point the config to them.

## Generating publication figures
1) Copy and adjust the config:
```
cp graph_config.json graph_config.local.json   # or edit graph_config.json directly
```
Fill in the CSV paths for ROC/AUC, long-range ROC, runtime/memory, beam sweep, etc. Example defaults point to `result/bpRNA_multifasta/ROC_analysis/*.csv`, `result/time_memory/...`, and `result/sarscov2/ensemble_energy/default/beam_sweep.csv`.

2) Run the figure pipeline:
```
uv run python -m scripts.generate_publication_figs --config graph_config.json
```
Outputs land in `graphs/`, including:
- `accuracy_auc*.{pdf,png}`: CapR vs LinCapR AUC by context (main + long-range + Turner1999 variants).
- `runtime_memory/*.pdf`: runtime & peak memory scatterplots (CapR, LinCapR, optionally CapR-cubic if enabled).
- `beam_sweep*.{pdf,png}`: SARS-CoV-2 beam sweep (time/memory and ensemble ΔG separated).
- Dataset summary and multiloop-unpaired figures when inputs are provided.

### Runtime/Memory with CapR-cubic
- Measurement script (span ≈ unlimited):
```
uv run python -m scripts.analysis.time_memory.measure_capr_cubic \
  --fasta data/processed/bprna/bpRNA_multifasta.fasta \
  --max-span 2000 \
  --output result/time_memory/capr_cubic/bprna_capr_cubic.csv \
  --capr-bin CapR/CapR
```
  - Default bins: 100–600 (6 seqs each), 600–1000 (3 seqs each). Tweak via `--bins`.
  - On macOS, `gtime` is preferred; otherwise `/usr/bin/time -l` is used.
- Enable plotting: set `runtime_memory_include_capr_cubic: true` and add the cubic CSV directory to `runtime_memory_input_dirs` in `graph_config.json`.

### Long-range stem AUC
Long-range ROC/AUC expects `roc_curves_longrange.csv` and `macro_auc_longrange.csv` with `context_range` values like `Stem_150+`. Positive counts per threshold are printed when plotting.

### Energy model subsets
Config fields `accuracy_allowed_energies`, `accuracy_turner1999_allowed_energies`, and their long-range counterparts let you restrict plots to specific energy models (e.g., only Turner1999 LinCapR).

## Notes
- `README.jp.md` contains the original Japanese guide.
- If a required CSV is missing, the corresponding figure is skipped with a log message.
- Colors: CapR (orange), LinCapR (sky blue), CapR-cubic (dark brown). CapR-cubic is shown without a regression line or beam legend entry.
