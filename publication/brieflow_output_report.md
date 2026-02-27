# Brieflow Output Assessment Report

**Output Directory:** `analysis/brieflow_output`

## 1. Preprocessing

- **SBS Tiles**: 13,488
- **Phenotype Tiles**: 20,760
- **Total Tiles**: 34,248

## 2. SBS (Sequencing by Synthesis)

- **Total Segmented Cells**: 18,092,081

### Barcode Mapping

- **Cells with 1 barcode**: 12,026,258 (66.5%)
- **Cells with 1+ barcodes**: 13,210,438 (73.0%)
- **Cells with 1 gene**: 6,928,656 (38.3%)
- **Cells with 1+ genes**: 7,257,370 (40.1%)

## 3. Phenotype

- **Total Cells**: 15,647,674
- **Feature Count**: 1,111

## 4. Merge & Aggregate Pipeline

### Inputs

- **SBS cells** (with 1+ barcodes): 13,210,438
- **Phenotype cells**: 15,647,674

### Spatial Matching

- **Candidate matches** (fast_merge): 13,954,151

### Deduplication

- **Initial matches**: 13,954,151 (5,873,374 with genes)
- **After SBS dedup** (1 phenotype per SBS cell): 13,279,390 (5,778,804 with genes)
- **After phenotype dedup** (1 SBS per phenotype cell): 12,113,480 (5,251,366 with genes)
- *Removed by deduplication*: 1,840,671

### Aggregate Filtering

- **Input to aggregate**: 12,113,480
- **After filtering** (requires gene + quality): 5,198,240
- *Removed* (no gene / quality): 6,915,240

**Overall Retention:** 5,198,240 / 13,210,438 = **39.3%**

### Aggregation Results

**Dataset:** CeCl-all_ChCo-Hoescht_COX4_AGP_ConA

- **Total Cells**: 5,198,240
- **Distinct Perturbations**: 21,732
- **Control Perturbations**: 1,624
- **Median Cells/Perturbation**: 224
- **Cell Range/Perturbation**: 6 - 1,517
- **PC Features**: 794

## 5. Clustering

### Full Clustering (all perturbations)

**Hoescht_COX4_AGP_ConA/all** (21,730 genes)

- **Resolution 1.0**: 48 clusters, CORUM 10 (N/A), KEGG 7 (11.2x)
- **Resolution 5.0**: 454 clusters, CORUM 24 (5.1x), KEGG 20 (2.4x)
- **Resolution 10.0**: 781 clusters, CORUM 28 (2.9x), KEGG 24 (1.9x)
- **Resolution 20.0**: 1256 clusters, CORUM 29 (2.4x), KEGG 31 (2.4x)
- **Resolution 40.0**: 1959 clusters, CORUM 35 (6.1x), KEGG 28 (2.2x)
- **Resolution 60.0**: 2553 clusters, CORUM 36 (12.3x), KEGG 26 (5.3x)
- **Resolution 80.0**: 3012 clusters, CORUM 40 (4.6x), KEGG 28 (2.1x)
- **Resolution 100.0**: 3345 clusters, CORUM 44 (5.7x), KEGG 30 (3.4x)

**Best:** Resolution 100.0 with 74 enriched clusters

### Filtered Clustering (bootstrap-significant perturbations)

**Hoescht_COX4_AGP_ConA/all** (3,715 genes)

*Filter: z-score threshold=0.3, FDR threshold=0.5, mode=both*

- **Resolution 10.0**: 320 clusters, STRING 10.0% precision, CORUM 15 (N/A), KEGG 8 (N/A)
- **Resolution 11.0**: 334 clusters, STRING 10.1% precision, CORUM 14 (N/A), KEGG 9 (N/A)
- **Resolution 12.0**: 363 clusters, STRING 10.5% precision, CORUM 15 (N/A), KEGG 8 (N/A)
- **Resolution 13.0**: 388 clusters, STRING 10.5% precision, CORUM 16 (N/A), KEGG 9 (10.0x)
- **Resolution 14.0**: 415 clusters, STRING 11.5% precision, CORUM 15 (N/A), KEGG 9 (N/A)
- **Resolution 15.0**: 427 clusters, STRING 10.4% precision, CORUM 16 (N/A), KEGG 9 (10.0x)

**Best:** Resolution 13.0 with 25 enriched clusters

#### MozzareLLM Analysis

*Model: claude-sonnet-4-5-20250929, Resolution: 10.0*

- **Total clusters analyzed**: 320
- **High confidence**: 27
- **Medium confidence**: 168
- **Low confidence**: 125

- **Established genes**: 776
- **Novel gene roles**: 1,050
- **Uncharacterized genes**: 758

**Top biological processes:**

- Cell cycle progression and mitotic regulation (6 clusters)
- Cytoskeletal organization and cell morphology regulation (5 clusters)
- Chromatin remodeling and transcriptional regulation (4 clusters)
- Endoplasmic reticulum protein processing and quality control (4 clusters)
