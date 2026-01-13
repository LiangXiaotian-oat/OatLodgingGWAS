# OatLodgingGWAS

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Language](https://img.shields.io/badge/language-Python%20%7C%20R-green.svg)

This repository contains the source code and analysis scripts for the research paper:

**"Multi-environment genome-wide association study identifies stable QTLs and candidate genes for lodging resistance-related traits in oat (_Avena sativa_ L.)"**

&gt; **Abstract / Overview:**  
&gt; This study utilizes multi-environment GWAS to dissect the genetic architecture of lodging resistance in oat. The repository provides a complete workflow for statistical analysis, including haplotype block identification, phenotype-genotype association mapping, and visualization of GWAS results.

---

## 📂 Repository Structure

The scripts are numbered and organized according to the analytical workflow presented in the manuscript, covering phenotype statistics, core trait screening, GWAS, and functional gene verification.

| Number | File Name | Description |
| :--- | :--- | :--- |
| **01** | `01_phenotype_stats.R` | <br>Calculates descriptive statistics (Mean, SD, CV, Range) for phenotypic data. |
| **02** | `02_correlation_visualization_PerformanceAnalytics.R` | <br>Analyzes phenotypic correlations across environments using `PerformanceAnalytics`. |
| **03** | `03_correlation_visualization_GGally.R` | <br>Visualizes pairwise correlation matrices and distributions separated by environment using `GGally`. |
| **04** | `04_xgboost_feature_importance.py` |<br>Calculates feature importance of agronomic traits contributing to LS using `XGBoost` for core trait screening. |
| **05** | `05_normality_checks_visualization.R` | <br>Plots frequency distribution histograms with fitted normal curves for normality checks. |
| **06** | `06_Heritability_BLUP_Calculation.R` | <br>Estimates broad-sense heritability ($H^2$) and BLUP values using Linear Mixed Models (LMM). |
| **07** | `07_gapit_gwas.R` | <br>Performs GWAS using the FarmCPU model in `GAPIT` package across environments and for BLUP values. |
| **08** | `08_Manhattan_QQ_Plots.R` | <br>Generates Manhattan and Q-Q plots for individual environments and BLUP values. |
| **09** | `09_SNP_genotype_visualization.py` | <br>Visualizes KASP genotyping results (scatter plots) and gene structures. |
| **10** | `10_haplotype_verification.py` | <br>1. Combined boxplots for single SNP effects across all traits. 2. Pyramiding effect analysis with regression lines (validating additive genetic effects). |

---

### 🛠️ Dependencies 

Please ensure the following dependencies are installed:

#### R Requirements:
- `GAPIT` (for GWAS)
- `lme4` (for BLUP/Heritability)
- `PerformanceAnalytics`, `GGally` (for Correlation)
- `ggplot2`, `CMplot` (for Visualization)
- `tidyverse`, `data.table` (Data manipulation)

#### Python Requirements:
- `pandas`, `numpy` (Data processing)
- `scikit-learn`, `xgboost` (Machine Learning)
- `matplotlib`, `seaborn` (Visualization)
- `scipy`, `statsmodels` (Statistical testing)

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
