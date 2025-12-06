# OatLodgingGWAS

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Language](https://img.shields.io/badge/language-Python%20%7C%20R-green.svg)

This repository contains the source code and analysis scripts for the research paper:

**"Multi-environment genome-wide association study identifies stable QTLs and candidate genes for lodging resistance-related traits in oat (*Avena sativa* L.)"**

## 📂 Repository Structure

* **`python_scripts/`**: Contains Python scripts for:
    * Haplotype block construction and identification.
    * Phenotype-genotype association analysis (ANOVA + Tukey's HSD).
    * Visualization of haplotype structures and genotype distributions.
* **`r_scripts/`**: Contains R scripts for:
    * Genome-wide association study (GWAS) visualization (Manhattan & QQ plots).
    * Statistical analysis and other visualizations used in the manuscript.
* **`data/`**: Sample datasets to demonstrate the usage of the scripts.

## 🚀 Getting Started

### Prerequisites

* **Python**: pandas, numpy, scipy, matplotlib, seaborn, statsmodels
* **R**: ggplot2, CMplot (or other packages you used)

### Usage Example (Haplotype Analysis)

To perform haplotype analysis using the provided Python script:

```bash
cd python_scripts
python haplotype_analysis.py --vcf ../data/example.vcf --pheno ../data/traits.csv --out results
