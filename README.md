# Whitney OPS Screen Analysis (v1.0.0)

Genome-wide fixed-cell optical pooled CRISPR-KO screen from Kirby, Di Bernardo et al. 2026 (in preparation), processed with [brieflow v1.4.6](https://github.com/cheeseman-lab/brieflow).

## Screen Overview

| Parameter | Value |
|-----------|-------|
| Cell line | HeLa-TetR-Cas9 (A7) |
| Genes targeted | 20,553 (genome-wide) |
| Total guides | 41,906 (4 per gene, dual-guide constructs) |
| Nontargeting controls | 200 |
| Plates | 2 |
| SBS cycles | 12 |
| Phenotype markers | Hoechst (DNA), COX4 (mitochondria), Phalloidin (actin/Golgi/PM), Concanavalin A (ER) |
| Microscope | Cephla Squid+ |

## Analysis

Raw data was processed end-to-end with [brieflow](https://github.com/cheeseman-lab/brieflow), covering preprocessing, SBS decoding, phenotype feature extraction, merge, aggregate, and clustering steps. Configuration notebooks and Slurm scripts are in `analysis/`.

## Data

Raw imaging data and brieflow outputs are stored at `/lab/ops_analysis_hdd/cheeseman/whitney-analysis/`. Source data used for figure generation is in `analysis/source_data/`.

## Citation

> *Citation will be added upon publication.*
