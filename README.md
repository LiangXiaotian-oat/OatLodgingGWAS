# OatLodgingGWAS

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Language](https://img.shields.io/badge/language-Python%20%7C%20R-green.svg)

This repository contains the source code and analysis scripts for the research paper:

**"Multi-environment genome-wide association study identifies stable QTLs and candidate genes for lodging resistance-related traits in oat (*Avena sativa* L.)"**

## 📂 Repository Structure

* **`python_scripts/`**: Contains Python scripts for:
    * Haplotype block construction and identification.
    * Phenotype-genotype association analysis (One-way ANOVA + Tukey's HSD).
    * Visualization of haplotype structures and genotype distributions.
* **`r_scripts/`**: Contains R scripts for:
    * Genome-wide association study (GWAS) visualization (Manhattan & QQ plots).
    * Statistical analysis (Correlation analysis, Heritability & BLUP calculation).

## 🛠️ Usage

Since raw datasets are not included in this repository due to size/privacy constraints, please prepare your input files as follows to run the scripts:

### 1. Python Scripts (Haplotype Analysis)
**Input Requirements:**
* **VCF File**: Standard VCF format containing SNPs for the target region.
* **Phenotype File**: CSV format. 
    * Column 1: `SampleID` (Must match VCF sample names).
    * Column 2+: Trait values.

**Command:**
```bash
cd python_scripts
python haplotype_analysis.py --vcf <your_data.vcf> --pheno <your_traits.csv> --out <output_dir>
