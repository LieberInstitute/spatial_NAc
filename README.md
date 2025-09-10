spatialNAc
================

<!-- README.md is generated from README.Rmd. Please edit that file -->


<img src="https://github.com/LieberInstitute/spatial_NAc/blob/main/project_overview.png?raw=true" width="1000px" align="left" />


## Overview

Welcome to the `spatialNAc` project! This project involves paired
snRNA-seq and SRT (10x Visium) data as well as several interactive
websites, all of which you are publicly accessible for you to browse and
download.

In this project we studied spatially resolved and single nucleus
transcriptomics data from the human Nucleus Accumbens (NAc) from
postmortem human brain samples. From 10 neurotypical controls we
generated spatially-resolved transcriptomics data using using [10x
Genomics
**Visium**](https://www.10xgenomics.com/products/spatial-gene-expression)
across the anterior, intermediate, and posterior NAc. We also generated
single nucleus RNA-seq (**snRNA-seq**) data using [10x Genomics
**Chromium**](https://www.10xgenomics.com/products/single-cell-gene-expression)

This project involves the GitHub repository
[LieberInstitute/NAc](https://github.com/LieberInstitute/spatialNAc)

If you tweet about this website, the data or the R package please use
the <code>\#spatialNAc</code> hashtag. You can find previous tweets that
way as shown
<a href="https://twitter.com/search?q=%23spatialDLPFC&src=typed_query">here</a>.

Thank you for your interest in our work!

## Study Design

**Study design to generate paired single nucleus RNA-sequencing
(snRNA-seq) and spatially-resolved transcriptomic data across NAc**.
Tissue blocks containing the NAc were dissected from 10 neurotypical
adult donors (6 male, 4 female). Paired single-nucleus RNA sequencing
(snRNA-seq) and spatially-resolved transcriptomics (SRT) data were
generated from adjacent tissue sections using 10x Genomics Chromium and
Visium platforms. To capture the complete NAc, tissue blocks were scored
to align with the width of the Visium capture array and spatial
profiling was performed across 2–5 capture arrays per donor (n = 38
total). snRNA-seq was performed on the same tissue blocks on PI sorted
(PI+) and neuron-enriched (PI+ NeuN+) nuclei. Downstream analyses
included (i) Identification of transcriptionally distinct cell clusters
from snRNA-seq, (ii) Mapping spatial domains in SRT data, (iii)
Integration of snRNA-seq and SRT data via cell-type deconvolution to
resolve spatially-localized populations, (iv) Inference of disease
relevant ligand-receptor (LR)–based cell–cell communication networks,
and (v) Cross-species drug-response mapping by integrating
transcriptional programs derived from rodent datasets with human SRT
data.

## Interactive Websites

All of these interactive websites are powered by open source software,
namely:

- 🔍 [`spatialLIBD`](https://doi.org/10.1186/s12864-022-08601-w)
- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

We provide the following interactive websites, organized by dataset with
software labeled by emojis:

- Visium 
  - 🔍 [spatial NAC](https://interactive.libd.org/spatial_NAC/)
- iSEE
  - 👀 [snRNA_NAC](https://interactive.libd.org/snRNA_NAC/)

## Data Availability  

### Zenodo Archive

There is a Zenodo archive for this repository available at [10.5281/zenodo.17089020](https://doi.org/10.5281/zenodo.17089020)

### Processed Data  

A public globus endpoint was created and contains objects too large for github. It is available at [https://app.globus.org/file-manager?origin_id=1d521cdd-4319-4719-acae-69d9d1ddc843&origin_path=%2F](https://app.globus.org/file-manager?origin_id=1d521cdd-4319-4719-acae-69d9d1ddc843&origin_path=%2F)  

### Raw Data

Sequencing data has been uploaded to GEO with accession numbers [GSE307586](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307586) for spatial data and [GSE307587](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307587) for single nucleus data

## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/spatialNAc/issues](https://github.com/LieberInstitute/spatialNAc/issues)
and refrain from emailing us. Thank you again for your interest in our
work!

## Citing our work

## Internal

- JHPCE locations:
  - `/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc`
