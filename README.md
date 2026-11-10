# Gene_network_GBM

Analysis script for immune-related gene interaction perturbation networks in glioblastoma (GBM). It builds a sample-specific edge perturbation matrix from GBM tumors against normal brain (GTEx), defines four subtypes (C1-C4) by consensus clustering, validates them in independent cohorts, and characterizes the subtypes.

## Citation

Liu Z, Xu Y, Wang Y, Weng S, Xu H, Ren Y, Guo C, Liu L, Zhang Z, Han X. Immune-related interaction perturbation networks unravel biological peculiars and clinical significance of glioblastoma. *iMeta* 2023;2(3):e127. doi: [10.1002/imt2.127](https://doi.org/10.1002/imt2.127)

The script matches the study described in that abstract: an interaction-perturbation framework on TCGA/CGGA GBM RNA-seq and GTEx normal brain, four subtypes validated in public datasets and an in-house cohort.

## What the script does

`Gene_network_GBM.R` is one long script (about 2,700 lines) whose sections are marked with `#...-----` comment headers, in this order:

1. Edge perturbation matrix: `lzq_genenetwork()` compares GBM RNA-seq (`Meta_RNAseq`) with GTEx expression on a gene interaction network (`net.rda`); `lzq_selectfeature()` keeps 4,000 edges.
2. Consensus clustering with `ConsensusClusterPlus` (PAM, 1,000 repetitions, K up to 10), with K = 4 taken for the subtypes C1-C4; consensus heatmap and survival curves.
3. Validation in a microarray cohort (`Meta_Array`) and an in-house ZZU cohort (`ZGBM`): subtype signature genes are extracted, samples are assigned with nearest template prediction (`CMScaller::ntp`), followed by survival analysis and heatmaps.
4. SubMap comparison of subtypes across cohorts, and subtype distribution across the three cohorts (alluvial and proportion plots).
5. Functional characterization: GO analysis per subtype, GSVA with hallmark gene sets, ESTIMATE immune and stromal scores and tumor purity, proliferation score, immune cells and checkpoints, APS and TIS scores, MHC expression, and cancer-immunity cycle steps.
6. Immunotherapy response prediction with SubMap against published cohorts (GSE35640, GSE91061, GSE100797, GSE126044, GSE135222 and one labelled NATH).
7. Genomic features: mutation comparison across subtypes, IDH mutation, G-CIMP, 1p/19q co-deletion, and TERT promoter mutation (also in the ZZU cohort).

## Dependencies

R with the following packages, as attached by the script:

`tidyverse`, `data.table`, `readxl`, `xlsx`, `magrittr`, `reshape2`, `survival`, `survminer`, `timeROC`, `pROC`, `caret`, `glmnet`, `ConsensusClusterPlus`, `CMScaller`, `clusterProfiler`, `org.Hs.eg.db`, `GSVA`, `GSEABase`, `msigdbr`, `limma`, `IOBR`, `MCPcounter`, `maftools`, `ComplexHeatmap`, `circlize`, `pheatmap`, `ggplot2`, `ggpubr`, `ggsci`, `ggalluvial`, `ggcor`, `ggkm`, `ggquickeda`, `compareGroups`, `RColorBrewer`, `export`

## Data and files not included

The script reads inputs that are not in this repository, so it cannot be run as-is:

- `Database/`: GBM RNA-seq expression and clinical data (`Meta_RNAseq*`), microarray cohort (`Meta_Array*`), ZZU in-house cohort (`ZGBM_expr`, `ZGBM_clin`), GTEx expression (`GTEx_expr`), the gene annotation (`GRCh38_v39.gtf`) and the interaction network (`net`).
- `Resource/lzq_genenetwork.R`: helper functions used by the script (`lzq_genenetwork`, `lzq_selectfeature`, `lzq_get_ntp_template`, `lzq_COC_normalize`).
- `Reference/`: gene sets for the immunotherapy-related analyses (for example `07.immune/CIC_Marker.rda`).
- Output folders `RData/`, `Figure/` and `Table/` with their sub-folders; the script writes into them and creates only a few itself, so create them first.

Run the script from the project root that contains these folders.
