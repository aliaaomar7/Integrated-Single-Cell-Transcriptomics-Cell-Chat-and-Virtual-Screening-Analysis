# Single-Cell-Transcriptomics-Cell-Chat-and-Virtual-Screening-Analysis
Integrated Single-Cell Transcriptomics, Cell-Cell Communication Inference, and Virtual Screening Prioritise CXCR4 and Natural Compound Candidates for Experimental Testing in Idiopathic Pulmonary Fibrosis

> **Status:** Preprint — all findings are in silico predictions requiring experimental validation.

---

## Table of Contents

- [Overview](#overview)
- [Key Findings](#key-findings)
- [Repository Structure](#repository-structure)
- [Methods Summary](#methods-summary)
  - [1. scRNA-seq Re-analysis](#1-scrna-seq-re-analysis)
  - [2. Cell-Cell Communication Inference](#2-cell-cell-communication-inference)
  - [3. Structure-Based Virtual Screening](#3-structure-based-virtual-screening)
  - [4. Network Pharmacology](#4-network-pharmacology)
- [Data Availability](#data-availability)
- [Dependencies](#dependencies)
- [Reproducing the Analysis](#reproducing-the-analysis)
- [Results at a Glance](#results-at-a-glance)
  - [Top Virtual Screening Hits](#top-virtual-screening-hits)
  - [Signalling Pathway Changes](#signalling-pathway-changes)
- [Proposed Experimental Validation](#proposed-experimental-validation)
- [Limitations](#limitations)
- [Citation](#citation)
- [References](#references)
- [License](#license)

---

## Overview

Idiopathic pulmonary fibrosis (IPF) is a progressive, fatal interstitial lung disease with a median survival of 3–5 years and only two FDA-approved drugs (pirfenidone, nintedanib) that slow — but do not reverse — disease progression.

This repository contains all code and resources for a three-layer computational study that:

1. **Re-analyses public scRNA-seq data** (GSE132771; n=3 IPF, n=3 normal donors) to characterise cell-type-level transcriptomic changes in IPF.
2. **Infers cell-cell communication rewiring** using CellChat v2 to identify hyperactivated signalling pathways in IPF.
3. **Performs structure-based virtual screening** of 105 bioactive natural compounds against the CXCR4 crystal structure (PDB: 3ODU) using AutoDock Vina, prioritising candidates that exceed the binding affinity of AMD3100 (Plerixafor).

The convergence of all three analytical layers on the **CXCL/CXCR4 axis** as a central signalling hub in IPF is the principal hypothesis-generating contribution of this work.

---

## Key Findings

| Layer | Finding |
|---|---|
| **scRNA-seq** | Fibroblasts expanded +65% (12.2% → 20.1%); epithelial cells declined −67% (16.1% → 5.3%) in IPF |
| **scRNA-seq** | CTHRC1⁺ pathological myofibroblasts appear exclusively in IPF (replicates Adams et al. 2020) |
| **scRNA-seq** | 189 upregulated / 134 downregulated genes in IPF fibroblasts; top GO terms: ECM organisation, apoptosis regulation, unfolded protein response |
| **CellChat v2** | 697 inferred interactions in IPF vs 502 in normal (+38.8%); total interaction strength +36.7% |
| **CellChat v2** | Top gained pathways: SPP1 (10.0), CCL (9.1), **CXCL (9.0)**, HGF (8.0), FGF (7.0) |
| **CellChat v2** | MIF → CD74/CXCR4 co-receptor engagement prominent; fibroblast autocrine self-loop dominant in IPF |
| **Virtual screening** | 10/101 compounds (9.9%) exceeded AMD3100 affinity (−9.471 kcal/mol) |
| **Virtual screening** | Top 4 prioritised hits: Coptisine (−10.560), Genistein (−9.973), Glycyrrhetinic acid (−9.876), Licochalcone A (−9.645) kcal/mol |
| **Network pharmacology** | CXCR4 hub module (n=8 proteins): 25 observed vs 5 expected edges; PPI enrichment p = 2.27 × 10⁻¹⁰ |

---

## Repository Structure

```
.
├── data/
│   └── README_data.md          # Instructions for downloading GSE132771 from GEO
├── scrnaseq/
│   ├── 01_qc_integration.R     # Seurat v5 QC, CCA integration, UMAP, clustering
│   ├── 02_annotation.R         # Cell type annotation, DotPlot, FeaturePlots
│   ├── 03_de_enrichment.R      # Fibroblast DE analysis, GO/KEGG enrichment (enrichR)
│   └── figures/                # Output figures (Figure 1, Figure 2)
├── cellchat/
│   ├── 04_cellchat_v2.R        # CellChat v2 pipeline, pathway analysis, bubble plots
│   └── figures/                # Output figures (Figure 3)
├── docking/
│   ├── compounds/              # SDF/PDBQT files for 105 natural compound library
│   ├── receptor/               # Prepared CXCR4 receptor (3ODU, polar H + Gasteiger charges)
│   ├── 05_prepare_receptor.sh  # PyMOL + OpenBabel receptor preparation
│   ├── 06_run_vina.sh          # AutoDock Vina batch docking script
│   ├── 07_parse_results.py     # Score extraction, Lipinski filtering, tier classification
│   └── figures/                # Output figures (Figure 4)
├── network_pharmacology/
│   ├── 08_string_ppi.R         # STRING v12.0 PPI network construction and enrichment
│   └── figures/                # Output figures (Figure 6)
├── results/
│   ├── docking_scores_all101.csv
│   ├── top_hits_table1.csv
│   └── cellchat_pathway_summary.csv
├── environment/
│   ├── r_packages.R            # R package installation script
│   └── conda_env.yml           # Python/AutoDock Vina conda environment
└── README.md
```

---

## Methods Summary

### 1. scRNA-seq Re-analysis

- **Data:** GEO accession [GSE132771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771) (Adams et al. 2020, *Science Advances*)
- **Samples:** GSM3891627, GSM3891629, GSM3891631 (IPF); GSM3891621, GSM3891623, GSM3891625 (Normal)
- **Pipeline:** Seurat v5 | min.cells=3, min.features=200 | Per-sample QC (nFeature 500–4,000; nCount <20,000; MT% <10%) | Downsampled to 3,000 cells/donor (set.seed=42) → **17,847 cells post-QC**
- **Integration:** CCA with 2,000 variable features (VST); top 20 PCs; UMAP + clustering (resolution=0.1) → 14 clusters → 8 cell types
- **Annotation markers:**

  | Cell Type | Markers |
  |---|---|
  | Monocyte-derived Macrophages | S100A8, FCN1 |
  | Endothelial | CLDN5, PECAM1, EMCN |
  | Epithelial (AT2-dominant) | EPCAM, SFTPC |
  | NK Cells | KLRD1, PRF1 |
  | Fibroblasts | COL1A2, LUM, PDGFRA/B |
  | Dendritic Cells | CD1C, FCER1A |
  | Plasma B Cells | IGHG1, IGKC |
  | Mast Cells | TPSB2 |

- **DE analysis:** `FindMarkers()` — Wilcoxon rank-sum, |log2FC| >0.5, adj. p <0.05 (Benjamini-Hochberg)
- **Enrichment:** enrichR — GO Biological Process 2023, KEGG 2021 Human

### 2. Cell-Cell Communication Inference

- **Tool:** CellChat v2 (Jin et al. 2021, 2025)
- **Database:** CellChatDB.human
- **Key parameters:** `computeCommunProb(type='triMean', raw.use=TRUE)`, `filterCommunication(min.cells=10)`, permutation test (10,000 shuffles)
- **NMF patterns:** CON k=4; IPF k=5

### 3. Structure-Based Virtual Screening

- **Target:** CXCR4 crystal structure — PDB [3ODU](https://www.rcsb.org/structure/3ODU) (2.5 Å)
- **Receptor preparation:** PyMOL (remove water + co-crystal ligand IT1t) → OpenBabel v3.1 (polar H + Gasteiger charges, pH 7.4)
- **Compound library:** 105 bioactive natural compounds from PubChem (MMFF94-optimised 3D coordinates); 4 excluded (2 cytotoxic, 1 CNS-liable, 1 duplicate) → **101 compounds docked**
- **Docking:** AutoDock Vina 1.2.5 | Grid: 30×30×30 Å centred on IT1t pose (X=20.610, Y=−7.972, Z=71.068 Å) | exhaustiveness=32, 9 binding modes/compound
- **Controls:** AMD3100/Plerixafor (positive; CID 65015) | Pirfenidone (negative; CID 40632)
- **Drug-likeness:** Lipinski Rule-of-Five assessed for all hits

### 4. Network Pharmacology

- **Tool:** STRING v12.0 (confidence ≥ 0.700)
- **Hub proteins:** CD44, CXCR4, MIF, SPP1, CXCL12, CD74, ITGB1, CCL2
- **Hub scores:** Sum of z-scores of degree, betweenness, and closeness centrality

---

## Data Availability

| Resource | Source | Accession / Link |
|---|---|---|
| scRNA-seq raw data | NCBI GEO | [GSE132771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771) |
| CXCR4 crystal structure | RCSB Protein Data Bank | [PDB: 3ODU](https://www.rcsb.org/structure/3ODU) |
| Natural compound structures | PubChem | CIDs listed in `results/top_hits_table1.csv` |
| STRING PPI | STRING v12.0 | [string-db.org](https://string-db.org) |

Raw scRNA-seq data are not redistributed here — download directly from GEO using the accession above.

---

## Dependencies

### R (≥ 4.3)

```r
# Install via r_packages.R
packages <- c(
  "Seurat",        # v5
  "CellChat",      # v2
  "enrichR",
  "ggplot2",
  "dplyr",
  "patchwork"
)
```

### Python / Docking Environment

```bash
conda env create -f environment/conda_env.yml
conda activate ipf_docking
# Includes: autodock-vina 1.2.5, openbabel 3.1, rdkit, pandas, matplotlib
```

### External Tools

- [PyMOL](https://pymol.org) — receptor preparation
- [AutoDock Vina 1.2.5](https://vina.scripps.edu)
- [OpenBabel 3.1](https://openbabel.org)

---

## Reproducing the Analysis

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>

# 2. Download GEO data (requires GEOquery or manual download)
Rscript data/download_geo.R   # downloads GSE132771 CellRanger outputs

# 3. scRNA-seq QC, integration, annotation
Rscript scrnaseq/01_qc_integration.R
Rscript scrnaseq/02_annotation.R
Rscript scrnaseq/03_de_enrichment.R

# 4. CellChat v2 analysis
Rscript cellchat/04_cellchat_v2.R

# 5. Receptor preparation (requires PyMOL + OpenBabel)
bash docking/05_prepare_receptor.sh

# 6. Batch docking
conda activate ipf_docking
bash docking/06_run_vina.sh        # ~2–4 hours on 8-core CPU

# 7. Parse and rank results
python docking/07_parse_results.py

# 8. Network pharmacology
Rscript network_pharmacology/08_string_ppi.R
```

> **Tip:** Steps 1–4 require ~16 GB RAM. Step 6 (docking) is parallelisable; set `--cpu` flag in `run_vina.sh` to match available cores.

---

## Results at a Glance

### Top Virtual Screening Hits

| Compound | PubChem CID | Affinity (kcal/mol) | MW (Da) | cLogP | Lipinski RO5 | Pose Convergence |
|---|---|---|---|---|---|---|
| **Coptisine** | 72322 | **−10.560** | 320.2 | 2.81 | Pass | MEDIUM |
| **Genistein** | 65752 | **−9.973** | 270.2 | 3.10 | Pass | MEDIUM |
| **Glycyrrhetinic acid** | 265237 | **−9.876** | 470.4 | 3.35 | Pass | MEDIUM |
| **Licochalcone A** | 5316743 | **−9.645** | 338.4 | 4.12 | Pass | **HIGH** |
| AMD3100 *(positive ctrl)* | 65015 | −9.471 | 502.5 | 0.42 | Borderline | — |
| Pirfenidone *(negative ctrl)* | 40632 | −6.098 | 185.2 | 0.90 | Pass | — |

### Signalling Pathway Changes

**Gained in IPF:** SPP1, CCL, CXCL, HGF, FGF, THY1, PTN, MK, GAP, VISFATIN

**Lost in IPF:** IFN-II, IL10, OCLN, CADM, CALCR, OSM, IL16, DHEAS, Adrenaline, FLRT, ADGRL, PTPRM

---

## Proposed Experimental Validation

| Tier | Assay |
|---|---|
| **Tier 1 (essential)** | Radioligand displacement (IC₅₀) for top 3 compounds against CXCR4 |
| **Tier 1 (essential)** | IPF fibroblast proliferation assay ± TGF-β stimulation |
| **Tier 2 (supporting)** | β-arrestin recruitment assay |
| **Tier 2 (supporting)** | CXCL12-induced fibroblast migration assay |
| **Tier 2 (supporting)** | Collagen gel contraction assay |
| **Tier 3 (in vivo)** | Bleomycin-induced pulmonary fibrosis (mouse): hydroxyproline, Ashcroft score, qPCR for Col1a1/Acta2, BAL cell differentials |

---

## Limitations

- All findings are **in silico predictions** — no experimental validation has been performed.
- Small scRNA-seq cohort (n=3 per group); compositional differences should be interpreted descriptively.
- Rigid-receptor docking (AutoDock Vina) does not capture CXCR4 conformational flexibility; molecular dynamics simulations would strengthen binding mode predictions.
- ADMET profiling limited to Lipinski RO5; full ADMET assessment (e.g. SwissADME, pkCSM) is recommended before progressing candidates.
- Potential promiscuous binders or aggregate-based false positives can only be excluded experimentally.

---

## Citation

If you use this code or data, please cite:

```
[Author(s)]. Integrated Single-Cell Transcriptomics, Cell-Cell Communication Inference,
and Virtual Screening Prioritise CXCR4 and Natural Compound Candidates for Experimental
Testing in Idiopathic Pulmonary Fibrosis. Preprint, 2025.
```

---

## References

Key references underpinning this study:

1. Adams TS et al. Single-cell RNA-seq reveals ectopic and aberrant lung-resident cell populations in IPF. *Sci Adv.* 2020;6:eaba1983.
2. Jin S et al. Inference and analysis of cell-cell communication using CellChat. *Nat Commun.* 2021;12:1088.
3. Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell communication. *Nat Protoc.* 2025;20:180–222.
4. Eberhardt J et al. AutoDock Vina 1.2.0. *J Chem Inf Model.* 2021;61:3891–3898.
5. Szklarczyk D et al. The STRING database in 2023. *Nucleic Acids Res.* 2023;51:D638–D646.
6. Park SY et al. Genistein inhibits CXCR4-mediated migration. *Cell Biochem Funct.* 2014;32:647–655.
7. Li X et al. Licochalcone A inhibits TGF-β1-induced fibroblast-to-myofibroblast transition. *Phytomedicine.* 2022;95:153854.

Full reference list available in the preprint.

---

## License

This repository is released under the [MIT License](LICENSE). The scRNA-seq data (GSE132771) is subject to the original data use terms from NCBI GEO.
