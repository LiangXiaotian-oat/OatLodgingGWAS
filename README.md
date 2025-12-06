# OatLodgingGWAS

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Language](https://img.shields.io/badge/language-Python%20%7C%20R-green.svg)

This repository contains the source code and analysis scripts for the research paper:

**"Multi-environment genome-wide association study identifies stable QTLs and candidate genes for lodging resistance-related traits in oat (_Avena sativa_ L.)"**

&gt; **Abstract / Overview:**  
&gt; This study utilizes multi-environment GWAS to dissect the genetic architecture of lodging resistance in oat. The repository provides a complete workflow for statistical analysis, including haplotype block identification, phenotype-genotype association mapping, and visualization of GWAS results.

---

## 📂 Code Organization

All analysis scripts are located in the root directory of this repository. They are categorized by their analytical function:

### 🧬 Haplotype & Association Analysis (Python)
* **`haplotype_analysis.py`**: The core pipeline for haplotype analysis. It performs:
    * **Haplotype Construction**: Identifies haplotype blocks from VCF files.
    * **Statistical Testing**: Combinatorial analysis of SNPs and One-way ANOVA with Tukey's HSD post-hoc test.
    * **Visualization**: Generates haplotype structure heatmaps, genotype distribution plots, and phenotype boxplots.
* **`xgboost_feature_analysis.py`**: Machine learning-based feature ranking to evaluate the importance of phenotypic traits.

### 📊 GWAS & Statistical Visualization (R)
* **`gapit_gwas_demo.R`**: A standardized workflow for running GWAS using FarmCPU and BLINK models via the `GAPIT` package.
* **`manhattan_qq_plot_demo.R`**: Generates publication-quality Manhattan and QQ plots using `CMplot`.
* **`correlation_analysis_demo.R`**: Performs correlation matrix analysis and visualization using `GGally`.
* **`heritability_blup_analysis.R`**: Calculates Broad-Sense Heritability ($H^2$) and BLUP values using Linear Mixed Models (`lme4`).
* **`normality_check_demo.R`**: Checks data normality and plots distribution histograms.

#### R Scripts (Visualization & Stats)
The R scripts are designed to be modular. You can run them in RStudio or via the command line (`Rscript script_name.R`).

| Module | Required Input Columns |
|--------|------------------------|
| **GWAS Plots** | `SNP`, `Chromosome`, `Position`, `P.value` |
| **Correlation/Normality** | Matrix: rows = samples, cols = phenotypic traits |
| **Heritability** | `Line`, `Env`, `Rep`, `Trait` |

---

## 📝 Data Availability

The raw genotype (GBS) and phenotype data used in this study are available from the corresponding author upon reasonable request, or refer to the accession codes provided in the manuscript (if applicable).  
For demonstration purposes, the scripts include **synthetic data generation modules**, allowing users to test the pipeline functionality without needing external datasets.

---

## 🖊️ Citation

If you use this code, logic, or methodology in your research, please cite our paper:

&gt; Liang Xiaotian, et al. (2025).  
&gt; *Multi-environment genome-wide association study identifies stable QTLs and candidate genes for lodging resistance-related traits in oat (Avena sativa L.)*.  
&gt; [Journal Name]. DOI: [Insert DOI here]

---

## 📜 License

This project is open-sourced under the **MIT License**.  
See the [LICENSE](LICENSE) file for details.

## 📧 Contact

For any questions regarding the code, please contact Liang Xiaotian ;494382219@qq.com&gt;.

---

## 🛠️ Getting Started & Usage

Since raw datasets are not included in this repository due to size and privacy constraints, please prepare your input files according to the requirements below to run the scripts.

### 1. Prerequisites

Ensure you have the following dependencies installed:

- **Python ≥ 3.8**  
  `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `statsmodels`, `scikit-learn`, `xgboost`

- **R ≥ 4.0**  
  `GAPIT`, `CMplot`, `ggplot2`, `lme4`, `emmeans`, `GGally`, `PerformanceAnalytics`

### 2. Python Scripts (Haplotype Analysis)

The Python pipeline performs end-to-end haplotype analysis.

#### Input Data Requirements

- **VCF file**: standard VCF containing SNPs for the target gene region, e.g. `target_region.vcf`  
- **Phenotype file**: CSV format  
  - Column 1: `SampleID` (must match sample names in the VCF header)  
  - Column 2+: trait values (numeric)

#### Running the Analysis

```bash
# Navigate to the python directory
cd python_scripts

# Run the main analysis script
python haplotype_analysis.py \
  --vcf ../data/your_data.vcf \
  --pheno ../data/traits.csv \
  --out ../results_dir
