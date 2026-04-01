# 🧬 16S rRNA Metagenomic Data Analysis Pipeline

> **Thesis Chapter 3 — Method Development**
> Full DADA2-based amplicon sequencing pipeline for processing, analysing, and visualising 16S rRNA gene sequencing data across 96 samples from three experimental datasets.

---

## 📖 Overview

This R Markdown script (`16S_metagenomic_analysis.Rmd`) implements an end-to-end 16S rRNA amplicon sequencing analysis pipeline for **96 samples** spanning three method development experiments:

| Dataset | Samples | Objective |
|---------|---------|-----------|
| **Enrichment Testing** | 20 | Evaluate microbial dynamics during selective enrichment for foodborne pathogens |
| **Host-DNA Depletion Testing** | 60 | Assess methods for reducing chicken host DNA contamination |
| **DNA Extraction Kit Comparison** | 14 | Compare commercial DNA extraction kits for food microbiome samples |

The pipeline covers everything from raw FASTQ files through to publication-quality visualisations, including quality filtering, denoising, chimera removal, taxonomic classification, phyloseq object construction, alpha/beta diversity analysis, and targeted pathogen tracking.

---

## ⚙️ Prerequisites

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| [R](https://www.r-project.org/) | ≥ 4.2.0 | Core language |
| [RStudio](https://posit.co/download/rstudio-desktop/) | Any recent | Recommended IDE |

### R Packages

The script includes installation chunks. To install all required packages manually:

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("dada2", "phyloseq", "ShortRead", "ComplexHeatmap"))

# CRAN packages
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "stringr",
  "RColorBrewer", "scales", "reshape2",
  "patchwork", "cowplot", "aplot",
  "pheatmap", "circlize", "tibble", "rmarkdown"
))
```

---

## 🗂️ Required Input Files

### 1. Raw FASTQ Sequences

Paired-end Illumina FASTQ files (`.fastq.gz`) should be stored in a directory and the path updated in the `Set Path` chunk:

```r
# In chunk: Set Path
path <- "/path/to/your/fastq/files"   # ← CHANGE THIS
```

Files must follow the naming convention:
- Forward reads: `*_R1*.fastq.gz`
- Reverse reads: `*_R2*.fastq.gz`

Sample names are extracted from the forward read filenames automatically.

### 2. SILVA Reference Database

Taxonomic classification requires the **SILVA nr99 v138.2** training set. Update the path in the `Assign Taxonomy` chunk:

```r
# In chunk: Assign Taxonomy to species level
training <- "/path/to/silva_nr99_v138.2_toSpecies_trainset.fa"   # ← CHANGE THIS
```

The SILVA training set can be downloaded from the [DADA2 reference databases page](https://benjjneb.github.io/dada2/training.html).

---

## 🚀 How to Run

### Option 1 — Knit in RStudio (Recommended)

1. Open `16S_metagenomic_analysis.Rmd` in RStudio
2. Update the FASTQ file path and SILVA database path (see above)
3. Click **Knit** → **Knit to HTML Notebook**

### Option 2 — Run from the R Console

```r
rmarkdown::render("16S_metagenomic_analysis.Rmd")
```

> ⚠️ **Note:** The DADA2 error learning and denoising steps are computationally intensive. Running on a high-performance computing (HPC) cluster is recommended for large datasets. The script was developed on the NeSI HPC platform.

---

## 🔬 Pipeline Workflow

```
Raw FASTQ files
      │
      ▼
Quality Visualisation (plotQualityProfile)
      │
      ▼
Trimming & Filtering
  truncLen = c(230, 200) | maxEE = c(2,2) | maxN = 0 | truncQ = 2
      │
      ▼
Dereplication
      │
      ▼
Error Learning & DADA2 Denoising (forward + reverse)
      │
      ▼
Paired-end Merging
      │
      ▼
Sequence Table Construction
      │
      ▼
Chimera Removal (removeBimeraDenovo)
      │
      ▼
Read Tracking (input → filtered → denoisedF → denoisedR → merged → nonchim)
      │
      ▼
Taxonomic Classification (SILVA nr99 v138.2, Kingdom → Species)
      │
      ▼
Phyloseq Object Construction (ps)
      │
      ▼
Dataset-Specific Subsetting & Analysis
  ├── Dataset 1: Enrichment Testing (20 samples)
  ├── Dataset 2: Host-DNA Depletion (60 samples)
  └── Dataset 3: DNA Extraction Kits (14 samples)
```

---

## 📊 Key R Objects Produced

### Global (All 96 Samples)

| Object | Description |
|--------|-------------|
| `ps` | Master phyloseq object — all 96 samples |
| `seqtab.nochim` | ASV table after chimera removal |
| `taxa` | Taxonomic classification table (Kingdom → Species) |
| `track` | Read counts at each DADA2 pipeline step |
| `ps_species_level_all96` | All samples agglomerated at Genus level |
| `ps_top20_rel_abund_all96` | Top 20 most abundant genera, relative abundance |

### Dataset 1 — Enrichment Testing

| Object | Description |
|--------|-------------|
| `ps_enrichment_testing` | 20 enrichment samples with Broth + Incubation metadata |
| `ps_family_level_enrichment` | Family-level agglomeration |
| `ps_top20_families_enrichment` | Top 20 families |
| `ps_rvs_only` | RVS broth subset (5 samples) |
| `ps_rvs_bpw` | RVS + BPW broths only (10 samples) |
| `ps_all_broths` | All enrichment samples with cleaned metadata |

### Dataset 2 — Host-DNA Depletion

| Object | Description |
|--------|-------------|
| `ps_hostDNAdep` | 60 host-DNA depletion samples with full metadata |
| `ps_hostDNAdep.top20_family` | Top 20 families, normalised |
| `ps_no_cheesecloth` | Filtered + Filtered Centrifuge samples only |
| `ps_cheesecloth` | Cheesecloth samples only |
| `ps_cheesecloth_incubation` | Cheesecloth + enriched samples only |

### Dataset 3 — DNA Extraction Kits

| Object | Description |
|--------|-------------|
| `DNAkits2` | 14 DNA kit comparison samples |
| `DNAkits2.top20` | Top 20 taxa, relative abundance |
| `DNAkits2.top20_family` | Top 20 families |

---

## 📁 Output Files

All figures are saved to the **working directory** by default. To change the output location, update the `file=` or `filename=` argument in each `ggsave()` / `jpeg()` / `png()` call.

| File | Description |
|------|-------------|
| `taxa20_all96samples.png` | Top 20 genera across all 96 samples |
| `taxa_top20_all_samples_with-other.png` | Top 20 genera with "Other" group (all 96 samples) |
| `read_counts_changes_through_DADA2.png` | Heatmap of read retention through pipeline |
| `alpha_diversity_enrichment.png` | Shannon + Simpson diversity over incubation time |
| `plot_top_genera_all_broths.png` | Top 20 genera faceted by enrichment broth |
| `plot_top_family_all_broths.png` | Top 20 families faceted by enrichment broth |
| `plot_Campylobacter_BB+CAT.png` | *Campylobacter* vs other bacteria (all broths) |
| `plot_pathogen_enrichment_all_broths.png` | *Salmonella* + *E. coli* vs other bacteria |
| `plot_Campylobacter_BB_CAT_RA_(all-bacteria).png` | *Campylobacter* species in Bolton broths |
| `16S_Campylobacter_species_vs_other_allbacteria_Bolton.png` | *Campylobacter* species vs all bacteria (Bolton) |
| `Campylobacter_species_read_counts_Bolton_broths.csv` | Read count summary for *Campylobacter* species |
| `heatmap_genera_annotated.png` | ComplexHeatmap — relative abundance + read count overlaid |
| `heatmap_genera_sidebyside.png` | ComplexHeatmap — relative abundance + read counts side by side |
| `heatmap_genera_log_annotated.png` | ComplexHeatmap — log-scaled with read count text |
| `genera_abundance_summary.csv` | Summary table of genus-level abundance + read counts |
| `taxa20_hostDNAdepletion-wocheesecloth.jpg` | Host-DNA depletion: filtered samples |
| `taxa20_hostDNAdepletion-nocheesecloth_ra.jpg` | Host-DNA depletion: filtered samples, relative abundance |
| `taxa20_hostDNAdepletion-cheesecloth.jpg` | Host-DNA depletion: cheesecloth samples |
| `taxa20_hostDNAdepletion-justcheesecloth_ra.jpg` | Host-DNA depletion: cheesecloth, relative abundance |
| `taxa20_hostDNAdepletion-enriched.jpg` | Host-DNA depletion: enriched samples |
| `taxa20_hostDNAdepletion-enriched_ra.jpg` | Host-DNA depletion: enriched samples, relative abundance |
| `taxa20_DNAkits.jpg` | DNA kit comparison — top 20 families |
| `taxa20_DNAkits_ra.jpg` | DNA kit comparison — relative abundance |

---

## 🧪 Experimental Design Details

### Dataset 1: Enrichment Testing

Four selective enrichment broths were tested at five incubation timepoints (0, 4, 8, 12, 24 hours):

| Sample Code | Broth | Purpose |
|-------------|-------|---------|
| `BFBB` | BF Bolton Base | *Campylobacter*-selective base medium |
| `CAT` | BF Bolton Base + CAT | Cefoperazone-Amphotericin-Tetrathionate supplement |
| `BPW` | Buffered Peptone Water | Non-selective pre-enrichment |
| `RVS` | Rappaport-Vassiliadis Soya | *Salmonella*-selective |

**Target pathogens tracked:** *Campylobacter* (Campylobacterota) and *Salmonella* / *E. coli* (Enterobacteriaceae)

### Dataset 2: Host-DNA Depletion Testing

A 2 × 3 × 3 factorial design testing:
- **Depletion treatment:** Host-DNA depletion (saponin + DNase) vs No depletion
- **Filtration method:** Cheesecloth / 0.45 µm Filtered / Filtered + Centrifuge
- **Enrichment:** BFBB w CAT / RVS / No Incubation

Includes biological triplicates (1, 2, 3) and procedural blanks, with diluted and undiluted (ND) variants.

### Dataset 3: DNA Extraction Kit Comparison

Three commercial kits compared across chicken rinse, spiked, and unspiked (control) samples:

| Sample Code | Kit |
|-------------|-----|
| `PrMigDNA` / `PRMgDNA` | Presto Mini gDNA Bacteria Kit (Geneaid) |
| `QDnBT` | DNeasy Blood & Tissue Kit (Qiagen) |
| `ZymoMiPr` / `ZymoMIPr` | Quick-DNA Miniprep Kit (Zymo Research) |

---

## 🎨 Colour Schemes

| Element | Colour | Hex |
|---------|--------|-----|
| *Campylobacter* | Teal green | `#1b9e77` |
| *Salmonella* | Orange | `#d95f02` |
| *Escherichia coli* | Teal green | `#1b9e77` |
| "Other" taxa | Light grey | `grey85` |
| "Unassigned" taxa | Mid grey | `grey65` |
| Top taxa | ColorBrewer `Paired` palette | — |
| Read depth heatmap | White → Dark Blue | `#23238c` |

---

## ⚠️ Known Issues & Notes

- **`combined_plot_cowplot` chunk** references `plot_abundance_families` and `plot_heatmap_reads` — these objects are not defined in the script and the chunk will error if run. They appear to be placeholders from an earlier version.
- **`aplot` combination chunk** references `pheatmap` as an object, which will error — `pheatmap()` returns invisibly and cannot be piped into `insert_right()` this way.
- **`ps_enrichment_testing.top20_enrichment_testing`** is referenced in several late chunks but is not defined in the script — these chunks will error if run without first creating this object.
- **`df_long`** is used without being defined in some chunks — ensure the long-format dataframe is created before running those sections.
- **File paths** are hardcoded to a NeSI HPC directory (`/nesi/nobackup/massey04083/...`). These must be updated before running on a different system.
- The script uses `multithread=FALSE` for `filterAndTrim()` for Windows compatibility. If running on Linux/macOS or HPC, set `multithread=TRUE` to speed up processing.

---

## 📋 Script Structure at a Glance

```
16S_metagenomic_analysis.Rmd
│
├── Package Installation & Loading
├── DADA2 Pipeline
│   ├── Set file paths
│   ├── Quality visualisation
│   ├── Trimming & filtering
│   ├── Dereplication
│   ├── Error learning & denoising
│   ├── Paired-end merging
│   ├── Sequence table & chimera removal
│   ├── Read tracking
│   └── Taxonomic assignment (SILVA)
│
├── Phyloseq Object Construction
│   ├── Alpha diversity (all 96 samples)
│   ├── Bray-Curtis NMDS (all 96 samples)
│   └── Top 20 genera plots (all 96 samples)
│
├── Read Depth Heatmap (DADA2 pipeline steps)
│
├── Dataset 1: Enrichment Testing Analysis
│   ├── Sample subsetting & metadata creation
│   ├── Alpha diversity (Shannon + Simpson)
│   ├── Bray-Curtis NMDS
│   ├── Top 20 genera — all broths (faceted)
│   ├── Top 20 families — all broths (faceted)
│   ├── ComplexHeatmap (relative abundance + read counts)
│   ├── Campylobacter vs other bacteria
│   ├── Salmonella + E. coli vs other bacteria
│   ├── Campylobacter species — Bolton broths
│   └── Paired statistical analysis (0h vs 24h)
│
├── Dataset 2: Host-DNA Depletion Analysis
│   ├── Sample subsetting & metadata creation
│   ├── Alpha diversity
│   ├── Filtered/Centrifuge samples — bar plots
│   ├── Cheesecloth samples — bar plots
│   └── Enriched cheesecloth samples — bar plots
│
└── Dataset 3: DNA Extraction Kit Analysis
    ├── Sample subsetting & metadata creation
    ├── Top 20 families — bar plots
    └── Relative abundance plots
```

---

## 📝 Citation

If you use or adapt this pipeline, please cite the key tools:

- **DADA2:** Callahan et al. (2016) *Nature Methods* 13:581–583. [doi:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- **phyloseq:** McMurdie & Holmes (2013) *PLOS ONE* 8(4):e61217. [doi:10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)
- **SILVA database:** Quast et al. (2013) *Nucleic Acids Research* 41:D590–D596. [doi:10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)

---

## 👩‍🔬 Author

Developed as part of a PhD thesis.
**Chapter 3 — Method Development: 16S rRNA Metabarcode Sequencing**

---

*Last updated: 2026*
