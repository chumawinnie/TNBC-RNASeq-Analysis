# TNBC-RNASeq-Analysis 🔬🧬📊💻🧪📈✨
Molecular subtyping of Triple Negative Breast Cancer (TNBC) using RNA-Seq analysis and immune deconvolution to identify novel biomarkers for early diagnosis and effective prognosis.
# TNBC-RNASeq-Analysis

## Project Overview

This project focuses on the molecular subtyping of Triple Negative Breast Cancer (TNBC) using RNA-Seq analysis and immune deconvolution. The aim is to identify novel biomarkers for early diagnosis and effective prognosis of TNBC, which is characterized by the absence of the three most common receptors targeted in breast cancer treatment.

---

## Objectives

1. Perform RNA-Seq analysis to identify differentially expressed genes (DEGs) in TNBC and Estrogen Receptor Positive (ERP) breast cancer samples.
2. Conduct functional enrichment analysis for Gene Ontology (GO) terms and KEGG pathways.
3. Analyze immune cell infiltration using deconvolution methods to understand the immune microenvironment in TNBC.

---

## Data

### Source

The RNA-Seq data used for this project was obtained from the E-MTAB-4993 dataset, containing:

- **63 samples**:
  - ERP: 51 samples
  - TNBC: 12 samples

### Preprocessing Tools

- **FASTQC**: For quality control of raw sequencing data.
- **FASTP**: For trimming and filtering low-quality reads.
- **HISAT2**: For alignment to the reference genome.
- **SamBamba**: For deduplication.
- **StringTie**: For transcript assembly and quantification.

---

## Methodology

1. **Quality Control and Preprocessing**:

   - Quality checks and trimming of raw reads.
   - Alignment to the reference genome (GRCh38).

2. **Differential Expression Analysis**:

   - Performed using tools such as StringTie and Ballgown.
   - Criteria for DEGs: p-value < 0.05, logFC > 1.5 or < -1.5.

3. **Functional Enrichment**:

   - **Gene Ontology (GO)**: Biological Process (BP), Cellular Component (CC), Molecular Function (MF).
   - **KEGG Pathways**: Identification of pathways enriched in upregulated and downregulated DEGs.

4. **Immune Cell Infiltration**:

   - Deconvolution using quantiseq and MCPCounter methods.

---

## Key Findings

1. **Differential Expression Analysis**:

   - Identified 15,947 DEGs initially, filtered to 9,771 DEGs with annotations.
   - Top 10 upregulated genes: FOXA1, RHOB, AR, CMBL, AGR2, ESR1, TFF3, SYBU, CBLC, DNALI1.
   - Top 10 downregulated genes: CENPW, EN1, A2ML1, TMSB15A, FOXC1, KRT16, SLC7A5, CDK6, MELTF, CA9.

2. **Gene Ontology Analysis**:

   - Significant enrichment in biological processes related to immune response and cancer biology.

3. **KEGG Pathways Analysis**:

   - Pathways associated with cell proliferation, apoptosis, and immune regulation were enriched.

4. **Immune Cell Infiltration**:

   - Deconvolution revealed distinct immune profiles in TNBC compared to ERP samples.

---

## Repository Structure

```
TNBC-RNASeq-Analysis/
├── scripts/                # Custom Python and R scripts
├── data/                   # Links to raw data sources
├── results/                # Output files (figures, tables, etc.)
└── README.md               # Project documentation
```

---

## How to Use

### Prerequisites

- **Python**: Version 3.8 or later
- **R**: Version 4.0 or later
- Required libraries and tools:
  - Python: `pandas`, `matplotlib`, `seaborn`
  - R: `DESeq2`, `Ballgown`, `EnrichR`
  - Command-line tools: `FASTQC`, `HISAT2`, `SamBamba`, `StringTie`

### Steps

1. Clone the repository:

   ```bash
   git clone https://github.com/chumawinnie/TNBC-RNASeq-Analysis.git
   ```

2. Navigate to the repository:

   ```bash
   cd TNBC-RNASeq-Analysis
   ```

3. Run the preprocessing pipeline:

   ```bash
   python scripts/fastqc.py
   python scripts/hisat2_SE.py
   ```

4. Perform differential expression and enrichment analysis:

   ```bash
   Rscript scripts/Deseq2-kallisto.r
   Rscript scripts/enrichment-analysis.r
   ```

5. Generate immune infiltration profiles:

   ```bash
   Rscript scripts/sleuth.r
   ```

---

## Authors

- **Chuma Winnie Obiora**

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.

---

## Acknowledgments

- The authors acknowledge the use of publicly available datasets and tools for RNA-Seq analysis.
- Special thanks to collaborators and mentors who contributed to the project.

