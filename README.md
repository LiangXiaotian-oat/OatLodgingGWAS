# OatLodgingGWAS

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Language](https://img.shields.io/badge/language-Python%20%7C%20R-green.svg)

This repository contains the source code and analysis scripts for the research paper:

**"Multi-environment genome-wide association study identifies stable QTLs and candidate genes for lodging resistance-related traits in oat (_Avena sativa_ L.)"**

&gt; **Abstract / Overview:**  
&gt; This study utilizes multi-environment GWAS to dissect the genetic architecture of lodging resistance in oat. The repository provides a complete workflow for statistical analysis, including haplotype block identification, phenotype-genotype association mapping, and visualization of GWAS results.

---

## 📂 Repository Structure / 代码库结构

本项目的代码按照论文的分析流程进行编号和组织，涵盖了从表型统计、核心性状筛选、GWAS 分析到候选基因功能验证的全过程。

The scripts are numbered and organized according to the analytical workflow presented in the manuscript, covering phenotype statistics, core trait screening, GWAS, and functional gene verification.

| 编号 | 文件名 (File Name) | 描述 (Description) |
| :--- | :--- | :--- |
| **01** | `01_phenotype_stats.R` | 计算表型数据的描述性统计量（均值、标准差、CV、极值等）。<br>Calculates descriptive statistics (Mean, SD, CV, Range) for phenotypic data. |
| **02** | `02_correlation_visualization_PerformanceAnalytics.R` | 使用 `PerformanceAnalytics` 包进行多环境表型相关性分析。<br>Analyzes phenotypic correlations across environments using `PerformanceAnalytics`. |
| **03** | `03_correlation_visualization_GGally.R` | 使用 `GGally` 包绘制分环境的性状相关性矩阵和分布图。<br>Visualizes pairwise correlation matrices and distributions separated by environment using `GGally`. |
| **04** | `04_xgboost_feature_importance.py` | 基于 Python `XGBoost` 模型计算各农艺性状对倒伏评分 (LS) 的特征重要性，用于筛选核心性状。<br>Calculates feature importance of agronomic traits contributing to LS using `XGBoost` for core trait screening. |
| **05** | `05_normality_checks_visualization.R` | 绘制表型数据的频次分布直方图并拟合正态曲线，进行正态性检验。<br>Plots frequency distribution histograms with fitted normal curves for normality checks. |
| **06** | `06_Heritability_BLUP_Calculation.R` | 基于混合线性模型 (LMM) 计算广义遗传力 ($H^2$) 和最佳线性无偏预测值 (BLUP)。<br>Estimates broad-sense heritability ($H^2$) and BLUP values using Linear Mixed Models (LMM). |
| **07** | `07_gapit_gwas.R` | 调用 `GAPIT` 包中的 FarmCPU 模型进行多环境及 BLUP 值的全基因组关联分析。<br>Performs GWAS using the FarmCPU model in `GAPIT` package across environments and for BLUP values. |
| **08** | `08_Manhattan_QQ_Plots.R` | 绘制单个环境及综合环境（BLUP）的曼哈顿图和 Q-Q 图。<br>Generates Manhattan and Q-Q plots for individual environments and BLUP values. |
| **09** | `09_SNP_genotype_visualization.py` | 可视化 KASP 标记的分型结果（散点图）及基因结构图。<br>Visualizes KASP genotyping results (scatter plots) and gene structures. |
| **10** | `10_haplotype_verification.py` | 绘制候选基因优异/非优异单倍型的表型差异箱线图（验证多效性）。<br>Plots boxplots showing phenotypic differences between haplotypes (verifying pleiotropic effects). |
| **11** | `11_pyramiding_effect_analysis.py` | 分析优异等位基因累加数量与表型之间的线性回归关系（聚合效应）。<br>Analyzes the linear regression between the number of superior alleles and phenotypes (pyramiding effect). |

---

### 🛠️ Dependencies / 依赖库

请确保安装了以下 R 包和 Python 库以运行上述代码：
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
