#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: 11_pyramiding_effect_analysis.py
Description:
    This script evaluates the cumulative effect (pyramiding effect) of superior alleles 
    across multiple candidate loci identified in GWAS.

    Scientific Goals:
    1. Identification: Determine the "Superior Allele" for each SNP based on trait direction 
       (e.g., lower is better for disease index, higher is better for yield).
    2. Quantification: Calculate the 'Pyramiding Score' (count of superior alleles) for each accession.
    3. Validation: Perform linear regression to verify the additive effect of stacking these alleles.

    Workflow:
    1. Data Integration: Merge VCF genotypes with multi-environment phenotype data.
    2. Directionality Check: Auto-detect superior alleles (Ref vs. Alt) per trait.
    3. Visualization:
       - Comparative Boxplots: Phenotypic distribution of haplotypes.
       - Single-Locus Regression: Additive effect of individual SNPs.
       - Pyramiding Regression: Linear relationship between Superior Allele Count (0-N) and Phenotype.
    4. Output: High-resolution plots (PNG/PDF) and detailed statistical tables.

Author:      Liang Xiaotian
Email:       494382219@qq.com
Date:        2025/12/20
License:     MIT
Dependencies: pandas, matplotlib, seaborn, scipy, numpy
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# ==========================================
#           User Configuration
# ==========================================

CONFIG = {
    # File Paths
    "VCF_FOLDER": ".",
    "PHENOTYPE_FILE": "phe386.csv",
    "OUTPUT_DIR": "11_Pyramiding_Analysis_Results",

    # Trait Definitions
    "TRAIT_GROUPS": ["LS", "TIL", "TID", "TIBR", "TIPS", "TIWT"],
    
    # Traits where LOWER values indicate better performance (Superior Allele logic)
    "LOWER_IS_BETTER_TRAITS": ["LS", "TIL"],

    # Environment/Location Prefixes
    "ENV_LIST": ["24WJ", "24BC", "25WJ", "25BC", "25CZ"],

    # Candidate SNPs for Pyramiding
    # Structure: "Label_on_Plot": {"file": "VCF_Name", "pos": [Position]}
    "TARGET_SNPS": {
        "S1C_329788453": {"file": "AVESA.00400a.r2.1Cg0002078", "pos": [329788453]},
        "S2D_5458489":   {"file": "AVESA.00400a.r2.2Dg0000076", "pos": [5458489]},
        "S2D_5674979":   {"file": "AVESA.00400a.r2.2Dg0000080", "pos": [5674979]},
        "S3D_44003226":  {"file": "AVESA.00400a.r2.3Dg0000747", "pos": [44003226]},
        "S3D_44003255":  {"file": "AVESA.00400a.r2.3Dg0000747", "pos": [44003255]}
    }
}

# ==========================================
#           Helper Functions
# ==========================================

def save_plot_dual_formats(filename_base, save_dir):
    """Saves the current matplotlib figure in both PNG (High-Res) and PDF (Vector) formats."""
    png_path = os.path.join(save_dir, f"{filename_base}.png")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    
    pdf_path = os.path.join(save_dir, f"{filename_base}.pdf")
    plt.savefig(pdf_path, bbox_inches='tight')
    print(f"   [SAVED PLOT] {filename_base}")

def load_phenotype(file_path):
    """Loads phenotype data, supporting both CSV and TSV formats."""
    print(f"[INFO] Loading phenotype file: {file_path}")
    try:
        df = pd.read_csv(file_path)
        if df.shape[1] < 2: 
            df = pd.read_csv(file_path, sep='\t')
    except Exception as e:
        print(f"[ERROR] Failed to load phenotype file: {e}")
        return None
    
    # Standardize first column to 'SampleID'
    df.rename(columns={df.columns[0]: 'SampleID'}, inplace=True)
    df['SampleID'] = df['SampleID'].astype(str)
    return df

def get_snp_data_from_vcfs(vcf_folder, target_snps_config):
    """Parses VCF files to extract genotypes for target loci."""
    all_snp_data = []
    
    # Optimize file access: Group targets by filename
    file_map = {}
    for label, info in target_snps_config.items():
        fname = info['file']
        if fname not in file_map: file_map[fname] = []
        file_map[fname].append((label, info['pos']))

    vcf_files = glob.glob(os.path.join(vcf_folder, "*.vcf"))
    
    for vcf_file in vcf_files:
        base_name = os.path.basename(vcf_file).replace('.vcf', '')
        if base_name not in file_map: continue
        
        target_configs = file_map[base_name]
        
        with open(vcf_file, 'r', encoding='utf-8', errors='ignore') as f:
            sample_names = []
            for line in f:
                line = line.strip()
                if line.startswith('##'): continue
                if line.startswith('#CHROM'):
                    parts = line.split('\t')
                    sample_names = parts[9:]
                    continue
                
                parts = line.split('\t')
                try: pos = int(parts[1])
                except ValueError: continue
                
                for label, pos_list in target_configs:
                    if pos in pos_list:
                        ref = parts[3]
                        alts = parts[4].split(',')
                        alleles_map = {0: ref}
                        for idx, alt in enumerate(alts):
                            alleles_map[idx + 1] = alt
                        
                        format_parts = parts[8].split(':')
                        try: gt_idx = format_parts.index('GT')
                        except ValueError: continue

                        for i, sample_data in enumerate(parts[9:]):
                            sample_id = sample_names[i]
                            gt_str = sample_data.split(':')[gt_idx]
                            gt_val = gt_str.replace('|', '/').split('/')
                            
                            bases = []
                            numeric_val = 0
                            valid = True
                            
                            # Parse Genotype (e.g., 0/0, 0/1)
                            for val in gt_val:
                                if val == '.': 
                                    bases.append('N')
                                    valid = False
                                else:
                                    try: 
                                        v = int(val)
                                        bases.append(alleles_map.get(v, 'N'))
                                        numeric_val += v
                                    except: 
                                        bases.append('N')
                                        valid = False
                            
                            if 'N' in bases: continue
                            
                            unique_bases = sorted(list(set(bases)))
                            hap = unique_bases[0] if len(unique_bases) == 1 else '/'.join(unique_bases)
                            
                            all_snp_data.append({
                                'SampleID': sample_id,
                                'SNP_Name': label,
                                'Haplotype': hap,
                                'Allele_Count': numeric_val if valid else np.nan
                            })
                            
    return pd.DataFrame(all_snp_data)

def reshape_phenotype(pheno_df, trait_groups, env_list):
    """Reshapes phenotype data from wide to long format."""
    melted_data = []
    value_cols = [c for c in pheno_df.columns if c != 'SampleID']
    
    for col in value_cols:
        matched_trait = None
        for trait in trait_groups:
            if col.endswith(trait):
                matched_trait = trait
                break
        
        if matched_trait:
            env = col.replace(matched_trait, '')
            if env in env_list:
                temp_df = pheno_df[['SampleID', col]].rename(columns={col: 'Value'})
                temp_df['Environment'] = env
                temp_df['Trait_Type'] = matched_trait
                melted_data.append(temp_df)
                
    if not melted_data: return pd.DataFrame()
    return pd.concat(melted_data, ignore_index=True)

# ==========================================
#           Scoring & Logic
# ==========================================

def calculate_superior_allele_count_single(snp_df, trait):
    """
    Standardizes single SNP allele count to 0 or 1.
    1 represents the Superior Allele based on trait direction.
    """
    group_0 = snp_df[snp_df['Allele_Count'] == 0]['Value']
    group_2 = snp_df[snp_df['Allele_Count'] == 2]['Value']
    
    if len(group_0) == 0 or len(group_2) == 0:
        return snp_df['Allele_Count'] / 2 # Fallback
    
    mean_0 = group_0.mean()
    mean_2 = group_2.mean()
    
    # Check direction preference
    is_lower_better = trait in CONFIG["LOWER_IS_BETTER_TRAITS"]
    
    if is_lower_better:
        alt_is_superior = mean_2 < mean_0
    else:
        alt_is_superior = mean_2 > mean_0
        
    if alt_is_superior:
        # Alt (2) is superior -> map to 1
        return snp_df['Allele_Count'] / 2
    else:
        # Ref (0) is superior -> map to 1
        return 1 - (snp_df['Allele_Count'] / 2)

def calculate_pyramiding_score_total(df_trait_env, trait):
    """
    Calculates the cumulative Pyramiding Score (0 to N) for each sample.
    Also extracts specific haplotype bases for the summary table.
    """
    snp_list = df_trait_env['SNP_Name'].unique()
    scores = []
    
    is_lower_better = trait in CONFIG["LOWER_IS_BETTER_TRAITS"]
    
    for snp in snp_list:
        snp_data = df_trait_env[df_trait_env['SNP_Name'] == snp]
        mean_0 = snp_data[snp_data['Allele_Count'] == 0]['Value'].mean()
        mean_2 = snp_data[snp_data['Allele_Count'] == 2]['Value'].mean()
        
        if pd.isna(mean_0) or pd.isna(mean_2): continue
        
        # Determine superior allele
        if is_lower_better:
            alt_is_superior = mean_2 < mean_0
        else:
            alt_is_superior = mean_2 > mean_0
            
        if alt_is_superior:
            score_map = {0: 0, 1: 0.5, 2: 1} # Alt is Superior
        else:
            score_map = {0: 1, 1: 0.5, 2: 0} # Ref is Superior
            
        for idx, row in snp_data.iterrows():
            if not pd.isna(row['Allele_Count']):
                s = score_map.get(row['Allele_Count'], 0)
                scores.append({'SampleID': row['SampleID'], 'Score': s})
                
    if not scores: return pd.DataFrame()
    score_df = pd.DataFrame(scores)
    
    # 1. Sum scores
    aggregated = score_df.groupby('SampleID')['Score'].sum().reset_index()
    aggregated.rename(columns={'Score': 'Total_Superior_Alleles'}, inplace=True)
    
    # 2. Get phenotypic values
    pheno_subset = df_trait_env[['SampleID', 'Value']].drop_duplicates()
    
    # 3. Pivot specific haplotypes for the output table
    haplotypes_wide = df_trait_env.pivot_table(
        index='SampleID', 
        columns='SNP_Name', 
        values='Haplotype', 
        aggfunc='first'
    ).reset_index()
    
    # 4. Merge
    temp_df = pd.merge(aggregated, pheno_subset, on='SampleID', how='inner')
    final_df = pd.merge(temp_df, haplotypes_wide, on='SampleID', how='left')
    
    return final_df

# ==========================================
#           Plotting Functions
# ==========================================

def plot_boxplot_only(df_plot, existing_snps, trait, env, save_dir, palette_dict):
    """Plot A: Comparative Boxplots with Significance."""
    n_snps = len(existing_snps)
    fig, axes = plt.subplots(1, n_snps, figsize=(2.0 * n_snps + 1, 5), sharey=True)
    if n_snps == 1: axes = [axes]
    plt.subplots_adjust(wspace=0)

    for i, snp in enumerate(existing_snps):
        ax = axes[i]
        snp_data = df_plot[df_plot['SNP_Name'] == snp].copy()
        
        # Sort based on trait direction
        mean_vals = snp_data.groupby('Haplotype')['Value'].mean()
        is_lower_better = trait in CONFIG["LOWER_IS_BETTER_TRAITS"]
        
        if is_lower_better:
            sorted_haps = mean_vals.sort_values(ascending=True).index.tolist()
        else:
            sorted_haps = mean_vals.sort_values(ascending=False).index.tolist()
            
        hap_counts = snp_data['Haplotype'].value_counts()
        new_labels = [f"{h}\n(n={hap_counts[h]})" for h in sorted_haps]
        
        sns.boxplot(x='Haplotype', y='Value', data=snp_data, order=sorted_haps, 
                    palette=palette_dict, showfliers=False, ax=ax, width=0.6)
        sns.stripplot(x='Haplotype', y='Value', data=snp_data, order=sorted_haps, 
                      color='black', alpha=0.3, size=3, ax=ax)
        
        # Add T-test stars
        if len(sorted_haps) >= 2:
            d1 = snp_data[snp_data['Haplotype'] == sorted_haps[0]]['Value']
            d2 = snp_data[snp_data['Haplotype'] == sorted_haps[1]]['Value']
            if len(d1) > 1 and len(d2) > 1:
                t, p = stats.ttest_ind(d1, d2, equal_var=False)
                if p < 0.001: star = "***"
                elif p < 0.01: star = "**"
                elif p < 0.05: star = "*"
                else: star = "ns"
                ax.text(0.5, 0.9, star, ha='center', va='center', 
                        transform=ax.transAxes, fontsize=12, fontweight='bold')

        ax.set_xticklabels(new_labels, rotation=0, fontsize=9)
        ax.set_xlabel(snp, fontsize=11, fontweight='bold', rotation=45)
        if i == 0: ax.set_ylabel(f"Phenotype ({trait})", fontsize=12)
        else:
            ax.set_ylabel("")
            ax.tick_params(left=False)
            ax.grid(False, axis='y')

    plt.suptitle(f"{trait} | {env} (Haplotype Comparison)", fontsize=14, y=0.98)
    
    # Legend
    legend_elements = [Patch(facecolor=color, edgecolor='black', label=hap) 
                       for hap, color in palette_dict.items()]
    axes[-1].legend(handles=legend_elements, title="Haplotype", 
                    bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    save_plot_dual_formats(f"{trait}_{env}_1_Boxplot", save_dir)
    plt.close()

def plot_regression_single(df_plot, existing_snps, trait, env, save_dir):
    """Plot B: Single Locus Additive Effect."""
    n_snps = len(existing_snps)
    fig, axes = plt.subplots(1, n_snps, figsize=(2.0 * n_snps + 1, 5), sharey=True)
    if n_snps == 1: axes = [axes]
    plt.subplots_adjust(wspace=0)

    for i, snp in enumerate(existing_snps):
        ax = axes[i]
        snp_data = df_plot[df_plot['SNP_Name'] == snp].copy()
        valid_reg = snp_data.dropna(subset=['Allele_Count', 'Value'])
        valid_reg['Superior_Count'] = calculate_superior_allele_count_single(valid_reg, trait)
        
        jitter = np.random.normal(0, 0.03, size=len(valid_reg))
        ax.scatter(valid_reg['Superior_Count'] + jitter, valid_reg['Value'], 
                   color='gray', alpha=0.5, s=15)
        
        if len(valid_reg['Superior_Count'].unique()) > 1:
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                valid_reg['Superior_Count'], valid_reg['Value'])
            x_vals = np.array([0, 1])
            y_vals = intercept + slope * x_vals
            ax.plot(x_vals, y_vals, color='red', linewidth=2, linestyle='--')
            
            # Stats (Top of Plot)
            reg_text = f"$R^2$={r_value**2:.2f}\n$P$={p_value:.1e}"
            ax.text(0.5, 1.02, reg_text, ha='center', va='bottom', transform=ax.transAxes, 
                    fontsize=9, color='darkred', fontweight='bold')
        else:
            ax.text(0.5, 0.5, "No Seg", ha='center', transform=ax.transAxes)

        ax.set_xlabel(snp, fontsize=11, fontweight='bold', rotation=45)
        ax.set_xticks([0, 1])
        ax.set_xlim(-0.3, 1.3)
        if i == 0: ax.set_ylabel(f"Phenotype ({trait})", fontsize=12)
        else:
            ax.set_ylabel("")
            ax.tick_params(left=False)
            
    plt.suptitle(f"{trait} | {env} (Single Locus Effect)", fontsize=14, y=1.1)
    save_plot_dual_formats(f"{trait}_{env}_2_Reg_Single", save_dir)
    plt.close()

def plot_regression_pyramiding(agg_df, trait, env, save_dir):
    """Plot C: Cumulative Pyramiding Effect."""
    df_plot = agg_df.dropna(subset=['Total_Superior_Alleles', 'Value'])
    if df_plot.empty: return

    plt.figure(figsize=(6, 6))
    ax = plt.gca()
    
    jitter = np.random.normal(0, 0.1, size=len(df_plot))
    ax.scatter(df_plot['Total_Superior_Alleles'] + jitter, df_plot['Value'], 
               color='#4c72b0', alpha=0.6, s=20, edgecolor='white', linewidth=0.5)
    
    if len(df_plot['Total_Superior_Alleles'].unique()) > 1:
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            df_plot['Total_Superior_Alleles'], df_plot['Value'])
        
        x_vals = np.array([df_plot['Total_Superior_Alleles'].min(), 
                           df_plot['Total_Superior_Alleles'].max()])
        y_vals = intercept + slope * x_vals
        ax.plot(x_vals, y_vals, color='#c44e52', linewidth=2.5, linestyle='-')
        
        # Stats (Outside Top)
        stats_text = f"$R^2$ = {r_value**2:.2f}, $P$ = {p_value:.1e}"
        ax.text(0.5, 1.02, stats_text, ha='center', va='bottom', transform=ax.transAxes,
                fontsize=12, fontweight='bold', color='black')
    
    ax.set_title(f"{trait} | {env}\nPyramiding Effect", fontsize=14, y=1.08)
    ax.set_xlabel("Number of Superior Alleles", fontsize=12)
    ax.set_ylabel(f"Phenotype ({trait})", fontsize=12)
    ax.set_xticks(range(0, 7)) # Ensure integers on X axis
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    save_plot_dual_formats(f"{trait}_{env}_3_Reg_Pyramiding", save_dir)
    plt.close()

# ==========================================
#           Main Execution
# ==========================================

def main():
    if not os.path.exists(CONFIG["OUTPUT_DIR"]): os.makedirs(CONFIG["OUTPUT_DIR"])
    
    print("[INFO] Step 1: Loading Data...")
    pheno_raw = load_phenotype(CONFIG["PHENOTYPE_FILE"])
    if pheno_raw is None: return
    
    pheno_long = reshape_phenotype(pheno_raw, CONFIG["TRAIT_GROUPS"], CONFIG["ENV_LIST"])
    
    print("[INFO] Step 2: Extracting SNP Genotypes...")
    geno_long = get_snp_data_from_vcfs(CONFIG["VCF_FOLDER"], CONFIG["TARGET_SNPS"])
    
    if geno_long.empty or pheno_long.empty:
        print("[ERROR] Data extraction failed. Please check your VCFs and Phenotype file.")
        return

    merged_df = pd.merge(geno_long, pheno_long, on='SampleID', how='inner')
    
    print("[INFO] Step 3: Starting Analysis Pipeline...")
    
    for trait in CONFIG["TRAIT_GROUPS"]:
        # Directory setup
        trait_dir = os.path.join(CONFIG["OUTPUT_DIR"], trait)
        tables_dir = os.path.join(trait_dir, "Tables")
        if not os.path.exists(trait_dir): os.makedirs(trait_dir)
        if not os.path.exists(tables_dir): os.makedirs(tables_dir)
        
        for env in CONFIG["ENV_LIST"]:
            subset = merged_df[(merged_df['Trait_Type'] == trait) & 
                               (merged_df['Environment'] == env)].copy()
            subset = subset.dropna(subset=['Value'])
            if subset.empty: continue
            
            # Ensure proper SNP ordering for plots
            snp_order = list(CONFIG["TARGET_SNPS"].keys())
            existing_snps = [s for s in snp_order if s in subset['SNP_Name'].unique()]
            if not existing_snps: continue
            
            print(f"   Processing: {trait} - {env}")

            # --- A. Save Single SNP Table ---
            subset[['SampleID', 'SNP_Name', 'Haplotype', 'Allele_Count', 'Value']]\
                .sort_values(['SNP_Name', 'SampleID'])\
                .to_csv(os.path.join(tables_dir, f"{trait}_{env}_SingleSNP_Data.csv"), index=False)
            
            # --- B. Plot 1: Boxplots ---
            all_haps = sorted(subset['Haplotype'].dropna().unique())
            base_palette = sns.color_palette("Set2", n_colors=len(all_haps))
            palette_dict = dict(zip(all_haps, base_palette))
            plot_boxplot_only(subset, existing_snps, trait, env, trait_dir, palette_dict)
            
            # --- C. Plot 2: Single Regression ---
            plot_regression_single(subset, existing_snps, trait, env, trait_dir)
            
            # --- D. Pyramiding Analysis ---
            agg_df = calculate_pyramiding_score_total(subset, trait)
            if not agg_df.empty:
                # Save Table with Haplotypes
                cols = ['SampleID', 'Total_Superior_Alleles', 'Value']
                snp_cols = [c for c in agg_df.columns if c not in cols]
                ordered_snp_cols = [s for s in snp_order if s in snp_cols]
                final_cols = cols + ordered_snp_cols
                
                agg_df[final_cols].sort_values('SampleID')\
                    .to_csv(os.path.join(tables_dir, f"{trait}_{env}_Pyramiding_Data.csv"), index=False)

                # Plot 3: Pyramiding Regression
                plot_regression_pyramiding(agg_df, trait, env, trait_dir)
            
    print(f"\n[SUCCESS] All results saved to: {CONFIG['OUTPUT_DIR']}")

if __name__ == '__main__':
    main()
