#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Haplotype Analysis Pipeline
---------------------------
Copyright (c) 2025 Liang Xiaotian <494382219@qq.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Description: 
    This script performs haplotype identification, combinatorial analysis, 
    and phenotype association mapping (ANOVA + Tukey's HSD) for specified SNP regions.
    It generates visualization plots including haplotype structures, genotype distributions, 
    and boxplots with significance markers.

Usage:
    python haplotype_analysis.py --vcf data.vcf --pheno traits.csv --out results_dir

Author: Liang Xiaotian
Email: 494382219@qq.com
Date: 2025/12/06
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import os
import itertools
import argparse
import sys
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# Configure Plotting Style
sns.set(style="ticks")
plt.rcParams['pdf.fonttype'] = 42  # For editable text in vector graphics
plt.rcParams['ps.fonttype'] = 42

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Haplotype Analysis and Phenotype Association Pipeline")
    
    parser.add_argument("--vcf", required=True, help="Path to the input VCF file (containing SNPs of interest)")
    parser.add_argument("--pheno", required=True, help="Path to the phenotype CSV file (Columns: SampleID, Trait1, Trait2...)")
    parser.add_argument("--out", default="Haplotype_Output", help="Output directory name (Default: Haplotype_Output)")
    parser.add_argument("--min_freq", type=float, default=0.05, help="Minimum haplotype frequency threshold (0.05 = 5%)")
    parser.add_argument("--max_combo", type=int, default=5, help="Maximum number of SNPs to combine for analysis (Default: 5)")
    parser.add_argument("--min_sample", type=int, default=3, help="Minimum sample size per group for statistical testing")
    
    return parser.parse_args()

def load_vcf_data(vcf_path, valid_samples):
    """
    Parses VCF file and returns a genotype matrix and SNP information.
    """
    print(f"[*] Loading VCF file: {vcf_path} ...")
    
    # Auto-detect header line
    header_line = 0
    try:
        with open(vcf_path, 'r', encoding='utf-8', errors='ignore') as f:
            for i, line in enumerate(f):
                if line.startswith("#CHROM"):
                    header_line = i
                    break
    except FileNotFoundError:
        sys.exit(f"Error: VCF file not found: {vcf_path}")
    
    try:
        df = pd.read_csv(vcf_path, sep='\t', skiprows=header_line)
    except Exception as e:
        sys.exit(f"Error reading VCF: {e}")

    if '#CHROM' in df.columns:
        df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)

    # Filter samples present in phenotype data
    vcf_samples = df.columns[9:]
    intersect_samples = [s for s in vcf_samples if s in valid_samples]
    
    print(f"    - Total SNPs detected: {len(df)}")
    print(f"    - Overlapping samples: {len(intersect_samples)}")
    
    if len(intersect_samples) == 0:
        sys.exit("Error: No overlapping samples found between VCF and Phenotype file.")

    snp_info = []
    snp_data = {}
    
    # Parse genotypes
    for _, row in df.iterrows():
        pos = str(row['POS'])
        ref = row['REF']
        alt = row['ALT']
        alts_list = alt.split(',')
        
        snp_info.append({'POS': pos, 'REF': ref, 'ALT': alt})
        
        col_data = {}
        for sample in intersect_samples:
            gt_str = str(row[sample]).split(':')[0]
            
            # Genotype calling logic
            base = 'N' # Default to Missing
            
            if '.' in gt_str or gt_str == 'nan':
                base = 'N'
            else:
                alleles = gt_str.replace('|', '/').split('/')
                if len(alleles) < 2:
                    base = 'N'
                elif alleles[0] != alleles[1]:
                    base = 'H' # Heterozygous
                elif alleles[0] == '0':
                    base = ref # Homozygous Ref
                else:
                    # Homozygous Alt (handle multi-allelic)
                    try:
                        alt_idx = int(alleles[0]) - 1
                        if 0 <= alt_idx < len(alts_list):
                            base = alts_list[alt_idx]
                        else:
                            base = 'N'
                    except:
                        base = 'N'
            col_data[sample] = base
            
        snp_data[pos] = col_data
        
    df_matrix = pd.DataFrame(snp_data).T # Rows: SNP, Cols: Sample
    df_info = pd.DataFrame(snp_info)
    
    return df_matrix.T, df_info # Return: Rows=Sample, Cols=SNP

def add_significance_bars(ax, data, x_col, trait, order):
    """
    Calculates Tukey's HSD and annotates significant differences on the plot.
    """
    try:
        tukey = pairwise_tukeyhsd(endog=data[trait], groups=data[x_col], alpha=0.05)
    except:
        return

    significant_pairs = []
    for row in tukey.summary().data[1:]:
        g1, g2, p_adj, reject = row[0], row[1], row[3], row[6]
        
        # Define significance labels
        if p_adj < 0.001: label = "***"
        elif p_adj < 0.01: label = "**"
        elif p_adj < 0.05: label = "*"
        else: label = "ns"
        
        significant_pairs.append(((g1, g2), label))

    if not significant_pairs: return

    # Plotting logic for bars
    y_max = data[trait].max()
    y_min = data[trait].min()
    y_range = y_max - y_min
    step = y_range * 0.12 
    current_height = y_max + step * 0.5
    
    group_map = {g: i for i, g in enumerate(order)}
    
    for (g1, g2), label in significant_pairs:
        x1 = group_map.get(g1)
        x2 = group_map.get(g2)
        if x1 is None or x2 is None: continue
        if x1 > x2: x1, x2 = x2, x1
        
        # Draw brackets
        bar_h = current_height
        bar_tips = bar_h - (y_range * 0.02)
        color = 'black' if label != "ns" else 'gray'
        lw = 1.0 if label != "ns" else 0.8
        
        ax.plot([x1, x1, x2, x2], [bar_tips, bar_h, bar_h, bar_tips], lw=lw, c=color)
        ax.text((x1 + x2) * 0.5, bar_h, label, ha='center', va='bottom', color=color, fontsize=9)
        
        current_height += step

    # Adjust Y-limit to fit bars
    ax.set_ylim(top=current_height + step)

def perform_anova_and_plot(data, x_col, trait, output_path, title_prefix, min_freq, min_n):
    """
    Performs ANOVA, generates boxplots with significance markers.
    Filters data based on missing values (N) and frequency.
    """
    # 1. Filter missing data
    clean_data = data[~data[x_col].astype(str).str.contains('N')].copy()
    clean_data = clean_data.dropna(subset=[trait])
    
    if clean_data.empty: return None

    # 2. Filter by frequency
    total_valid = len(clean_data)
    count_threshold = int(total_valid * min_freq)
    
    counts = clean_data[x_col].value_counts()
    valid_groups = counts[counts >= count_threshold].index.tolist()
    
    if not valid_groups: return None

    # 3. Sort groups by mean phenotype value (Descending)
    group_means = clean_data[clean_data[x_col].isin(valid_groups)].groupby(x_col)[trait].mean()
    valid_groups = sorted(valid_groups, key=lambda x: group_means.get(x, 0), reverse=True)
    
    clean_data = clean_data[clean_data[x_col].isin(valid_groups)]
    
    if len(valid_groups) < 2: return None

    # 4. ANOVA
    groups_val = [clean_data[clean_data[x_col] == g][trait].values for g in valid_groups]
    try:
        f, p = stats.f_oneway(*groups_val)
    except:
        return None

    sig_symbol = "ns"
    if p < 0.001: sig_symbol = "***"
    elif p < 0.01: sig_symbol = "**"
    elif p < 0.05: sig_symbol = "*"

    # 5. Plotting
    fig, ax = plt.subplots(figsize=(max(5, len(valid_groups)*1.5), 6 + len(valid_groups)*0.5))
    
    sns.boxplot(x=x_col, y=trait, data=clean_data, order=valid_groups, palette="Set2", showfliers=False, width=0.5, ax=ax)
    sns.stripplot(x=x_col, y=trait, data=clean_data, order=valid_groups, color='black', alpha=0.5, size=3, ax=ax)
    
    # Add significance bars
    add_significance_bars(ax, clean_data, x_col, trait, valid_groups)
    
    title_text = f"{title_prefix}\n(Freq >= {min_freq*100}%, N={total_valid})\nANOVA p={p:.2e} ({sig_symbol})"
    ax.set_title(title_text)
    ax.set_xlabel("Haplotype (Sorted by Mean)")
    ax.set_ylabel(trait)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    
    return {
        'Trait': trait,
        'Combination_Level': len(x_col.split('_')) if '_' in x_col else 1,
        'P_Value': p,
        'Significance': sig_symbol,
        'Sample_Size': len(clean_data),
        'Num_Haps': len(valid_groups)
    }

def plot_haplotype_structure(hap_map_df, snp_cols, output_file):
    """Visualizes the composition of haplotypes (Heatmap)."""
    if hap_map_df.empty: return
    
    matrix_data = []
    y_labels = []
    
    # Sort by count
    sorted_df = hap_map_df.sort_values('Count', ascending=False)
    
    for _, row in sorted_df.iterrows():
        label = f"{row['HapName']} (n={row['Count']})"
        y_labels.append(label)
        matrix_data.append(list(row['Sequence']))
        
    if not matrix_data: return

    hap_matrix = pd.DataFrame(matrix_data, index=y_labels, columns=snp_cols)
    
    # Base to Number mapping for coloring
    base_to_num = {'A': 1, 'T': 2, 'C': 3, 'G': 4, 'N': 0, 'H': 5, '-': 0}
    num_matrix = hap_matrix.applymap(lambda x: base_to_num.get(x, 0))
    
    plt.figure(figsize=(max(8, len(snp_cols)*0.8), max(5, len(y_labels)*0.6)))
    
    # Colors: N=Gray, A=Red, T=Blue, C=Green, G=Orange, H=Purple
    colors_map = {0:'#d9d9d9', 1:'#e41a1c', 2:'#377eb8', 3:'#4daf4a', 4:'#ff7f00', 5:'#984ea3'}
    present_vals = sorted(np.unique(num_matrix.values))
    present_colors = [colors_map[v] for v in present_vals]
    
    sns.heatmap(num_matrix, cmap=present_colors, 
                annot=hap_matrix, fmt='', annot_kws={"size": 10, "weight": "bold"},
                linewidths=1, linecolor='white', cbar=False)
    
    plt.title(f"Haplotype Structure ({len(snp_cols)} SNPs)", fontsize=14)
    plt.yticks(rotation=0)
    plt.xticks(rotation=45)
    
    # Custom Legend
    legend_labels = {0:'Miss', 1:'A', 2:'T', 3:'C', 4:'G', 5:'Het'}
    handles = [Patch(color=colors_map[v], label=legend_labels[v]) for v in present_vals]
    plt.legend(handles=handles, title="Base", bbox_to_anchor=(1.01, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_genotype_distribution(data, snp_cols, snp_info_df, output_file):
    """Plots the distribution of genotypes for each SNP (Grouped Bar Chart)."""
    print(f"    - Plotting Genotype Distribution: {output_file}")
    
    plot_data = []
    for snp in snp_cols:
        counts = data[snp].value_counts()
        try:
            ref_alt = snp_info_df[snp_info_df['POS'].astype(str) == str(snp)].iloc[0]
            label = f"{snp}\n({ref_alt['REF']}/{ref_alt['ALT']})"
        except:
            label = str(snp)
        
        for base in ['A', 'T', 'C', 'G', 'H', 'N']:
            count = counts.get(base, 0)
            if count > 0:
                plot_data.append({'SNP': label, 'Genotype': base, 'Count': count})
    
    if not plot_data: return

    df_plot = pd.DataFrame(plot_data)
    
    # Dynamic width
    fig_width = 5 if len(snp_cols) == 1 else max(8, len(snp_cols) * 1.5)
    plt.figure(figsize=(fig_width, 6))
    
    colors_dict = {'A':'#e41a1c', 'T':'#377eb8', 'C':'#4daf4a', 'G':'#ff7f00', 'H':'#984ea3', 'N':'#d9d9d9'}
    
    ax = sns.barplot(data=df_plot, x='SNP', y='Count', hue='Genotype', 
                     palette=colors_dict, edgecolor='black', linewidth=0.5)
    
    for container in ax.containers:
        ax.bar_label(container, padding=3, fontsize=9)

    plt.title("Genotype Distribution per SNP", fontsize=14)
    plt.xlabel("SNP Position (Ref/Alt)", fontsize=12)
    plt.ylabel("Sample Count", fontsize=12)
    plt.legend(title="Genotype", bbox_to_anchor=(1.02, 1), loc='upper left')
    
    if len(snp_cols) > 1: plt.xticks(rotation=45)
    else: plt.xticks(rotation=0)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

# ==========================================
# Main Execution Flow
# ==========================================

if __name__ == "__main__":
    args = parse_args()
    
    # Directories
    dir_single = os.path.join(args.out, "1_Single_SNP")
    dir_single_details = os.path.join(dir_single, "Detailed_Genotypes")
    dir_combined = os.path.join(args.out, "2_Combined_Haplotype")
    
    for d in [args.out, dir_single, dir_single_details, dir_combined]:
        if not os.path.exists(d):
            os.makedirs(d)
            
    # 1. Load Data
    try:
        pheno_df = pd.read_csv(args.pheno, dtype={0: str})
    except:
        pheno_df = pd.read_csv(args.pheno, sep='\t', dtype={0: str})
    
    pheno_df.rename(columns={pheno_df.columns[0]: 'SampleID'}, inplace=True)
    valid_samples = set(pheno_df['SampleID'])
    traits = pheno_df.columns[1:]
    
    # 2. Parse VCF
    df_matrix, df_snp_info = load_vcf_data(args.vcf, valid_samples)
    df_matrix.index.name = 'SampleID'
    df_matrix.reset_index(inplace=True)
    
    # 3. Merge
    full_data = pd.merge(df_matrix, pheno_df, on='SampleID')
    snp_cols = [c for c in df_matrix.columns if c != 'SampleID']
    num_snps = len(snp_cols)
    
    all_stats = []
    
    print(f"\n[*] Starting Analysis (Total SNPs: {num_snps})")
    print(f"    Filter: Freq >= {args.min_freq}, Min Sample >= {args.min_sample}")

    # ---------------------------------------------------------
    # Phase 1: Single SNP Analysis
    # ---------------------------------------------------------
    print("\n>>> Phase 1: Single SNP Analysis")
    
    # Distribution Plot
    plot_genotype_distribution(full_data, snp_cols, df_snp_info, 
                             os.path.join(dir_single, "Genotype_Distribution.png"))
    
    single_stats = []
    for snp in snp_cols:
        # Save details
        detail_df = full_data[['SampleID', snp]].rename(columns={snp: 'Genotype'})
        detail_df.to_csv(os.path.join(dir_single_details, f"SNP_{snp}_Samples.csv"), index=False)
        
        # Analysis
        ref_alt = df_snp_info[df_snp_info['POS']==snp].iloc[0]
        title = f"SNP: {snp} ({ref_alt['REF']}/{ref_alt['ALT']})"
        
        for trait in traits:
            out_path = os.path.join(dir_single, f"SNP_{snp}_{trait}.png")
            res = perform_anova_and_plot(full_data, snp, trait, out_path, title, args.min_freq, args.min_sample)
            if res:
                res['Level'] = 1
                res['SNP_Combination'] = snp
                all_stats.append(res)
                single_stats.append(res)
    
    if single_stats:
        pd.DataFrame(single_stats).to_csv(os.path.join(dir_single, "Single_SNP_Stats.csv"), index=False)

    # ---------------------------------------------------------
    # Phase 2: Combinatorial Analysis
    # ---------------------------------------------------------
    print("\n>>> Phase 2: Combinatorial Haplotype Analysis")
    
    limit_k = min(num_snps, args.max_combo)
    
    for k in range(2, limit_k + 1):
        level_name = f"Level_{k}_SNPs"
        if k == num_snps: level_name = "Level_Full_Haplotype"
        level_dir = os.path.join(args.out, level_name)
        if not os.path.exists(level_dir): os.makedirs(level_dir)
        
        print(f"    Processing Level {k} combinations...")
        
        for snp_subset in itertools.combinations(snp_cols, k):
            combo_name = "_".join([str(s) for s in snp_subset])
            if len(combo_name) > 50: combo_name = f"Combo_{k}SNPs_Hash{hash(combo_name) % 10000}"
            
            # Build Haplotype
            col_name = f'Hap_{combo_name}'
            full_data[col_name] = full_data[list(snp_subset)].apply(lambda x: ''.join(x), axis=1)
            
            # Haplotype naming & stats (filtering N for naming)
            valid_hap_rows = full_data[~full_data[col_name].str.contains('N')].copy()
            
            if valid_hap_rows.empty: continue
            
            hap_counts = valid_hap_rows[col_name].value_counts()
            hap_map = {seq: f"H{i+1:03d}" for i, (seq, _) in enumerate(hap_counts.items())}
            
            # Apply map to all data (including N)
            full_data['Hap_Label'] = full_data[col_name].map(hap_map).fillna('Rare/Missing')
            
            # Save Samples
            out_csv = os.path.join(level_dir, f"Samples_{combo_name}.csv")
            full_data[['SampleID', col_name, 'Hap_Label']].to_csv(out_csv, index=False)
            
            # Plot Structure (Only for Full or specific levels to save time)
            if k == num_snps or k == 2: 
                map_df = pd.DataFrame([{'HapName': v, 'Sequence': k, 'Count': hap_counts[k]} 
                                      for k, v in hap_map.items() if k in hap_counts])
                plot_haplotype_structure(map_df, list(snp_subset), os.path.join(level_dir, f"Structure_{combo_name}.png"))
            
            # Phenotype Analysis
            for trait in traits:
                title = f"{k} SNPs Combo\n{combo_name}"
                out_path = os.path.join(level_dir, f"{combo_name}_{trait}.png")
                res = perform_anova_and_plot(full_data, 'Hap_Label', trait, out_path, title, args.min_freq, args.min_sample)
                
                if res:
                    res['Level'] = k
                    res['SNP_Combination'] = combo_name
                    all_stats.append(res)

    # ---------------------------------------------------------
    # Summary
    # ---------------------------------------------------------
    if all_stats:
        df_res = pd.DataFrame(all_stats)
        cols = ['Level', 'Trait', 'P_Value', 'Significance', 'SNP_Combination', 'Num_Haps', 'Sample_Size']
        df_res = df_res[cols].sort_values(by=['Trait', 'P_Value'])
        
        df_res.to_csv(os.path.join(args.out, "Final_Summary_Stats.csv"), index=False)
        print(f"\n[Success] Analysis complete. Summary saved to {args.out}/Final_Summary_Stats.csv")
    else:
        print("\n[Info] Analysis complete but no significant associations found.")
