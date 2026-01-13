#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: 10_Figure6_Pleiotropy_Pyramiding.py
Description:
    This script generates the composite visualization for Figure 6 in the manuscript.
    It integrates genotypic (KASP) and phenotypic data to perform two key analyses:
    
    1. Single Locus Pleiotropy Validation (Left Panels):
       - Compares phenotypic distributions between alleles for candidate genes across multiple traits.
       - Performs statistical testing (t-test) to verify pleiotropic effects.
       
    2. Pyramiding Effect Analysis (Right Panels):
       - Quantifies the cumulative effect of superior alleles.
       - Performs linear regression to demonstrate additive genetic gains.

    Input Requirements:
    - Phenotype Data: Excel/CSV with SampleID and trait values.
    - Genotype Data: Excel/CSV with SampleID and SNP calls (e.g., A/T/C/G).

Author:      Liang Xiaotian
Email:       494382219@qq.com
Date:        2025/12/20
License:     MIT
Dependencies: pandas, matplotlib, seaborn, scipy, openpyxl
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# ==========================================
#            User Configuration
# ==========================================

CONFIG = {
    # Input File Paths (Ensure these exist in the working directory)
    "GENOTYPE_FILE": "table.xlsx",       # KASP genotyping results
    "PHENOTYPE_FILE": "phenomean.xlsx",  # Mean phenotype data
    
    # Output Directory
    "OUTPUT_DIR": "Analysis_Final_Results",
    
    # Target SNPs (Column names in Genotype File)
    "SNP_COLS": [
        "S1C_329788453", 
        "S2D_5458489", 
        "S2D_5674979", 
        "S3D_44003226", 
        "S3D_44003255"
    ],
    
    # Trait List (Column names in Phenotype File)
    "TRAIT_GROUPS": ["LS", "TIL", "TID", "TIBR", "TIPS", "TIWT"],
    
    # Logic Definition: Traits where LOWER values are better (e.g., Lodging Score)
    # Used to automatically determine which allele is "Superior"
    "LOWER_IS_BETTER_TRAITS": ["LS", "TIL"]
}

# ==========================================
#            Helper Functions
# ==========================================

def setup_plotting_style():
    """Configures matplotlib for publication-quality figures."""
    sns.set(style="ticks")
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    # Fallback fonts to avoid errors on servers without Arial
    plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif']

def save_plot_dual_format(filename, directory):
    """Saves plot in both PNG (preview) and PDF (vector) formats."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    png_path = os.path.join(directory, f"{filename}.png")
    pdf_path = os.path.join(directory, f"{filename}.pdf")
    
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"   [SAVED] {filename} (.png & .pdf)")

def is_lower_better(trait):
    """Checks if a lower value indicates a better phenotype for the given trait."""
    return trait in CONFIG["LOWER_IS_BETTER_TRAITS"]

def load_and_merge_data():
    """Loads genotype and phenotype files and merges them on SampleID."""
    print("[INFO] Loading data...")
    
    # Check if files exist
    if not os.path.exists(CONFIG["GENOTYPE_FILE"]) or not os.path.exists(CONFIG["PHENOTYPE_FILE"]):
        print(f"[ERROR] Input files not found. Please ensure '{CONFIG['GENOTYPE_FILE']}' and '{CONFIG['PHENOTYPE_FILE']}' are in the directory.")
        sys.exit(1)

    try:
        # distinct logic for xlsx vs csv could be added here, assuming xlsx per config
        df_geno = pd.read_excel(CONFIG["GENOTYPE_FILE"])
        df_pheno = pd.read_excel(CONFIG["PHENOTYPE_FILE"])
    except Exception as e:
        print(f"[ERROR] Failed to read files: {e}")
        sys.exit(1)

    # Clean column names and IDs
    df_geno.columns = df_geno.columns.str.strip()
    df_pheno.columns = df_pheno.columns.str.strip()
    df_geno['SampleID'] = df_geno['SampleID'].astype(str).str.strip()
    df_pheno['SampleID'] = df_pheno['SampleID'].astype(str).str.strip()

    # Merge
    merged = pd.merge(df_pheno, df_geno, on='SampleID', how='inner')
    print(f"[INFO] Data merged successfully. Total samples: {len(merged)}")
    return merged

# ==========================================
#            Analysis Functions
# ==========================================

def calculate_pyramiding_score(df, snp_cols, trait):
    """
    Calculates the 'Superior Count' for each sample based on trait directionality.
    Auto-detects superior alleles.
    """
    df_score = df[['SampleID', trait]].copy()
    df_score['Superior_Count'] = 0
    df_score = df_score.dropna(subset=[trait])
    
    lower_better = is_lower_better(trait)
    
    for snp in snp_cols:
        # Filter valid genotypes
        temp = df[['SampleID', snp, trait]].dropna()
        temp[snp] = temp[snp].astype(str).str.strip()
        temp = temp[temp[snp].isin(['A', 'T', 'C', 'G'])] # Only standard nucleotides
        
        if temp.empty: continue

        # Determine superior allele based on mean comparison
        means = temp.groupby(snp)[trait].mean()
        if len(means) < 2: continue 
        
        if lower_better:
            superior_allele = means.idxmin() # Min is best (e.g., LS)
        else:
            superior_allele = means.idxmax() # Max is best (e.g., TIBR)
            
        # Add score
        ids_with_superior = temp[temp[snp] == superior_allele]['SampleID']
        df_score.loc[df_score['SampleID'].isin(ids_with_superior), 'Superior_Count'] += 1
        
    return df_score

# ==========================================
#            Plotting Functions
# ==========================================

def plot_combined_single_snp(df, snp_cols, traits, output_dir):
    """
    Generates a grid of Boxplots: (Rows=Traits) x (Cols=SNPs).
    Validates pleiotropy by showing SNP effects across multiple traits.
    """
    print("[INFO] Generating Combined Single SNP Plot (Figure 6 Left)...")
    
    n_rows = len(traits)
    n_cols = len(snp_cols)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.5 * n_cols, 3.5 * n_rows), constrained_layout=True)
    
    for row_idx, trait in enumerate(traits):
        lower_better = is_lower_better(trait)
        
        for col_idx, snp in enumerate(snp_cols):
            ax = axes[row_idx, col_idx]
            
            # Prepare data
            data = df[['SampleID', snp, trait]].dropna()
            data.columns = ['SampleID', 'Genotype', 'Value']
            data['Genotype'] = data['Genotype'].astype(str).str.strip()
            data = data[data['Genotype'].isin(['A', 'T', 'C', 'G'])]
            
            if data.empty:
                ax.text(0.5, 0.5, "No Data", ha='center', transform=ax.transAxes)
                continue

            # Sort alleles: Non-superior first, Superior second (for visual consistency)
            geno_counts = data['Genotype'].value_counts()
            top_genos = geno_counts.index.tolist()[:2] # Take top 2 if >2 alleles exist
            
            means = data.groupby('Genotype')['Value'].mean()
            valid_means = means[means.index.isin(top_genos)]
            
            if lower_better:
                # Plot order: High Value (Bad) -> Low Value (Good)
                order = valid_means.sort_values(ascending=False).index.tolist()
            else:
                # Plot order: Low Value (Bad) -> High Value (Good)
                order = valid_means.sort_values(ascending=True).index.tolist()
            
            plot_data = data[data['Genotype'].isin(order)]
            labels = [f"{g}\n(n={geno_counts[g]})" for g in order]
            
            # Plotting
            palette = sns.color_palette("Set2", n_colors=len(order))
            sns.boxplot(x='Genotype', y='Value', data=plot_data, order=order, 
                        palette=palette, showfliers=False, ax=ax, width=0.5)
            sns.stripplot(x='Genotype', y='Value', data=plot_data, order=order, 
                          color='black', alpha=0.3, size=3, ax=ax)
            
            # T-test Significance
            if len(order) == 2:
                g1 = plot_data[plot_data['Genotype'] == order[0]]['Value']
                g2 = plot_data[plot_data['Genotype'] == order[1]]['Value']
                if len(g1) > 1 and len(g2) > 1:
                    t, p = stats.ttest_ind(g1, g2, equal_var=False)
                    
                    sig = "ns"
                    if p < 0.001: sig = "***"
                    elif p < 0.01: sig = "**"
                    elif p < 0.05: sig = "*"
                    
                    # Draw significance bar
                    y_max = plot_data['Value'].max()
                    y_min = plot_data['Value'].min()
                    rng = y_max - y_min if y_max != y_min else 1.0
                    h = rng * 0.1
                    y_line = y_max + h * 0.5
                    
                    ax.plot([0, 0, 1, 1], [y_line, y_line+h*0.2, y_line+h*0.2, y_line], lw=1, c='k')
                    ax.text(0.5, y_line+h*0.3, sig, ha='center', va='bottom', color='k', fontweight='bold')

            # Styling
            ax.set_xlabel("")
            ax.set_xticklabels(labels, fontsize=9)
            
            # Only show Y label on the first column
            if col_idx == 0:
                ax.set_ylabel(trait, fontsize=14, fontweight='bold')
            else:
                ax.set_ylabel("")
            
            # Only show SNP Title on the first row
            if row_idx == 0:
                ax.set_title(snp, fontsize=12, fontweight='bold')

    save_plot_dual_format("Figure6_Left_SingleSNP", output_dir)
    plt.close()

def plot_combined_pyramiding(df, snp_cols, traits, output_dir):
    """
    Generates a grid of Regression Plots: 2 Rows x 3 Cols.
    Shows the additive effect of accumulating superior alleles.
    """
    print("[INFO] Generating Pyramiding Effect Plot (Figure 6 Right)...")
    
    n_traits = len(traits)
    # Fixed layout for 6 traits
    n_rows, n_cols = 2, 3
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 12), constrained_layout=True)
    axes_flat = axes.flatten()
    
    for i, trait in enumerate(traits):
        if i >= len(axes_flat): break
        ax = axes_flat[i]
        
        # Calculate scores
        df_score = calculate_pyramiding_score(df, snp_cols, trait)
        
        # Prepare Plot Data
        plot_data = df_score[['Superior_Count', trait]].dropna()
        plot_data.columns = ['Count', 'Value']
        possible_counts = sorted(plot_data['Count'].unique())
        
        if len(possible_counts) == 0: continue

        # Boxplot
        sns.boxplot(x='Count', y='Value', data=plot_data, order=possible_counts, 
                    palette="viridis", showfliers=False, ax=ax, width=0.6)
        sns.stripplot(x='Count', y='Value', data=plot_data, order=possible_counts, 
                      color='black', alpha=0.3, size=4, ax=ax)
        
        # Linear Regression Line
        if len(possible_counts) > 1:
            slope, intercept, r_val, p_val, std_err = stats.linregress(plot_data['Count'], plot_data['Value'])
            
            # Plot line across the full range
            x_vals = np.array([min(possible_counts), max(possible_counts)])
            y_vals = intercept + slope * x_vals
            ax.plot(x_vals, y_vals, color='#c44e52', linewidth=2.5, linestyle='--', zorder=10)
            
            # Annotate Stats
            stats_text = f"$R^2$ = {r_val**2:.2f}\n$P$ = {p_val:.1e}"
            ax.text(0.5, 0.95, stats_text, ha='center', va='top', transform=ax.transAxes,
                    fontsize=12, fontweight='bold', 
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", alpha=0.8))

        # Styling
        ax.set_title(trait, fontsize=16, fontweight='bold')
        ax.set_xlabel("Superior Haplotypes Count", fontsize=10)
        ax.set_ylabel(f"Phenotype ({trait})", fontsize=10)
        
        # X-axis Count Labels
        counts_stats = plot_data['Count'].value_counts()
        new_labels = [f"{int(c)}\n(n={counts_stats.get(c,0)})" for c in possible_counts]
        ax.set_xticklabels(new_labels)

    # Hide unused subplots
    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].axis('off')
        
    save_plot_dual_format("Figure6_Right_Pyramiding", output_dir)
    plt.close()

# ==========================================
#            Main Execution
# ==========================================

def main():
    setup_plotting_style()
    
    # 1. Load Data
    df_merged = load_and_merge_data()
    
    # 2. Generate Left Panel (Single SNP Effects)
    plot_combined_single_snp(df_merged, CONFIG["SNP_COLUMNS"], CONFIG["TRAIT_GROUPS"], CONFIG["OUTPUT_DIR"])
    
    # 3. Generate Right Panel (Pyramiding Effects)
    plot_combined_pyramiding(df_merged, CONFIG["SNP_COLUMNS"], CONFIG["TRAIT_GROUPS"], CONFIG["OUTPUT_DIR"])
    
    print(f"\n[SUCCESS] Analysis complete. All figures saved to: '{CONFIG['OUTPUT_DIR']}'")

if __name__ == "__main__":
    main()
