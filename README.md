# OatLodgingGWAS

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Language](https://img.shields.io/badge/language-Python%20%7C%20R-green.svg)

This repository contains the source code and analysis scripts for the research paper:

**"Multi-environment genome-wide association study identifies stable QTLs and candidate genes for lodging resistance-related traits in oat (_Avena sativa_ L.)"**

&gt; **Abstract / Overview:**  
&gt; This study utilizes multi-environment GWAS to dissect the genetic architecture of lodging resistance in oat. The repository provides a complete workflow for statistical analysis, including haplotype block identification, phenotype-genotype association mapping, and visualization of GWAS results.

---

## 📂 Code Organization & Analysis Pipeline

All scripts are numbered according to the analysis workflow described in the manuscript. **Note:** All scripts include a built-in demo mode (`--demo` or synthetic data generation), allowing users to reproduce the pipeline without accessing the raw dataset.

### 📊 1. Phenotype Statistical Analysis
* **`01_phenotype_stats.R`**
    * Calculates descriptive statistics (Mean, SD, CV, Skewness, Kurtosis) for all traits across environments.
    * Pre-processes raw data by averaging replicates.
* **`02_correlation_visualization_PerformanceAnalytics.R`**
    * Performs Pearson correlation analysis within single environments.
    * Generates correlation matrices with significance levels and histograms using the `PerformanceAnalytics` package.
* **`03_correlation_visualization_GGally.R`**
    * Visualizes multi-environment trait correlations.
    * Uses `GGally` to create pairwise plots colored by environment to assess trait stability.
* **`05_normality_checks_visualization.R`**
    * Performs Shapiro-Wilk normality tests.
    * Generates high-resolution histograms with fitted normal distribution curves for quality control.

### 🤖 2. Machine Learning & Heritability
* **`04_xgboost_feature_importance.py`**
    * **Key Script:** Implements the XGBoost regression model described in Section 2.3.2.
    * Ranks phenotypic traits by their "Importance Score" (Gain) regarding Lodging Resistance.
    * Includes logic for grouping features (e.g., Morphological vs. Mechanical traits).
* **`06_Heritability_BLUP_Calculation.R`**
    * Fits Linear Mixed Models (LMM) using `lme4`.
    * Calculates Best Linear Unbiased Predictions (**BLUPs**) for GWAS input.
    * Estimates Broad-sense Heritability ($H^2$) using variance components.

### 🧬 3. Genome-Wide Association Study (GWAS)
* **`07_gapit_gwas.R`**
    * **Core Analysis:** Runs the GWAS pipeline using **GAPIT 3**.
    * Implements **FarmCPU** and **BLINK** models as specified in Section 2.4.
    * Handles synthetic HapMap genotype data generation for demonstration.
* **`08_Manhattan_QQ_Plots.R`**
    * Visualizes GWAS results using `CMplot`.
    * Generates publication-quality Rectangular and Circular Manhattan plots.

### 🔬 4. Candidate Gene & Haplotype Analysis
* **`09_SNP_genotype_visualization.py`**
    * Visualizes genotype frequency distributions for specific SNPs.
    * Generates stacked bar charts to inspect allele counts.
* **`10_haplotype_verification.py`**
    * **Validation Script:** Performs the haplotype analysis described in Section 2.5.2.
    * Parses VCF files to identify haplotype blocks.
    * Conducts **ANOVA** and **Tukey’s HSD test** to verify phenotypic differences between haplotypes.
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
