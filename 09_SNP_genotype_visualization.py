#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: 09_SNP_genotype_visualization.py
Description:
    This script visualizes SNP genotype distributions (Stacked/Grouped Bar Charts).
    It is designed to check the allele frequency or genotype counts for specific SNPs.
    
    Features:
    - Auto-generates synthetic data if input file is missing (Demo Mode).
    - Saves publication-quality figures (PNG & PDF) to standardized output dir.
    - Robust font handling for different operating systems.

Author:      Liang Xiaotian
Email:       494382219@qq.com
Date:        2025/12/09
License:     MIT
"""

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 1. Dependency Management
# ==========================================
def install_package(package):
    """Install a package using pip if import fails."""
    try:
        print(f"[Setup] Installing missing package: {package}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    except subprocess.CalledProcessError:
        print(f"[Error] Failed to install {package}. Please install manually.")

# Check required packages
required_packages = ['pandas', 'matplotlib', 'seaborn', 'openpyxl']
for lib in required_packages:
    try:
        __import__(lib)
    except ImportError:
        install_package(lib)

# ==========================================
# 2. Configuration
# ==========================================

# File Settings
# If this file exists, it will be loaded. Otherwise, synthetic data is generated.
INPUT_FILE = "single_snp_data.xlsx" 

# Output Settings (Standardized Project Structure)
OUTPUT_DIR = os.path.join("results", "09_SNP_Genotype_Vis")
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# Plot Settings
FIG_SIZE = (10, 6)
DPI_SETTING = 300

# Font Settings (Robust Fallback logic)
# Tries standard fonts to avoid errors on Linux/Servers
try:
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
except Exception:
    pass

# Color Settings (Standard Genotype Colors)
COLOR_PALETTE = {
    "AA": "#D6221E",  # Red
    "TT": "#3778AD",  # Blue
    "CC": "#4EA74A",  # Green
    "GG": "#EF7C1B",  # Orange
    "AT": "#9b59b6",  # Purple (Het)
    "NN": "gray"      # Missing
}

# ==========================================
# 3. Data Loading & Simulation logic
# ==========================================

def generate_synthetic_data():
    """Generates dummy SNP genotype data for demonstration."""
    print("[Data] Input file not found. Generating synthetic SNP data...")
    
    # Simulate 5 SNPs
    snps = [f"SNP_{i}" for i in range(1, 6)]
    genotypes = ["AA", "TT", "AT", "CC", "GG"]
    
    data = []
    for snp in snps:
        # Generate random counts for genotypes
        # Simulate different distributions for each SNP
        counts = np.random.randint(10, 100, size=len(genotypes))
        for gt, count in zip(genotypes, counts):
            data.append({"SNP": snp, "Genotype": gt, "Count": count})
            
    df = pd.DataFrame(data)
    
    # Save dummy data for user reference (Optional)
    dummy_path = os.path.join(OUTPUT_DIR, "demo_snp_data.csv")
    df.to_csv(dummy_path, index=False)
    print(f"[Data] Synthetic data saved to: {dummy_path}")
    
    return df

def load_dataset(file_path):
    """Loads data from file or generates synthetic data if missing."""
    if not os.path.exists(file_path):
        return generate_synthetic_data()

    try:
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        elif file_path.endswith(('.xls', '.xlsx')):
            df = pd.read_excel(file_path)
        else:
            print("[Error] Unsupported file format. Please use .csv or .xlsx")
            sys.exit(1)
        
        print(f"[Data] Successfully loaded: {file_path}")
        return df
    except Exception as e:
        print(f"[Error] Failed to read file: {e}")
        sys.exit(1)

# ==========================================
# 4. Plotting Function
# ==========================================

def plot_snp_distribution(df, output_folder):
    """Plots grouped bar chart for SNP genotypes."""
    
    # Validation: Ensure columns exist
    required_cols = {'SNP', 'Genotype', 'Count'}
    if not required_cols.issubset(df.columns):
        print(f"[Error] Dataframe missing required columns: {required_cols}")
        print(f"Found: {df.columns.tolist()}")
        return

    sns.set_style("ticks")
    plt.figure(figsize=FIG_SIZE)
    
    try:
        # Draw Barplot
        ax = sns.barplot(
            data=df,
            x="SNP",
            y="Count",
            hue="Genotype",
            palette=COLOR_PALETTE,
            edgecolor="black",
            linewidth=1.2
        )
    except Exception as e:
        print(f"[Error] Plotting failed: {e}")
        return

    # Styling
    plt.title("Genotype Distribution per SNP", fontsize=16, fontweight='bold', pad=20)
    plt.xlabel("SNP Marker", fontsize=14, fontweight='bold')
    plt.ylabel("Sample Count", fontsize=14, fontweight='bold')
    
    # Legend
    plt.legend(title="Genotype", title_fontsize='12', bbox_to_anchor=(1.02, 1), loc='upper left')
    
    # Add value labels on bars
    for container in ax.containers:
        ax.bar_label(container, padding=3, fontsize=10)

    plt.tight_layout()

    # Save outputs (Vector & Raster)
    png_path = os.path.join(output_folder, "09_SNP_Genotype_Plot.png")
    pdf_path = os.path.join(output_folder, "09_SNP_Genotype_Plot.pdf")
    
    plt.savefig(png_path, dpi=DPI_SETTING, bbox_inches='tight')
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    
    print(f"[Output] Plots saved to:\n  - {png_path}\n  - {pdf_path}")
    # plt.show() # Uncomment for interactive mode

# ==========================================
# 5. Main Execution
# ==========================================

if __name__ == "__main__":
    print("--- Starting SNP Visualization Script (09) ---")
    
    # 1. Load or Generate Data
    df = load_dataset(INPUT_FILE)
    
    # 2. Preview
    print("Data Preview (Head):")
    print(df.head())
    
    # 3. Plot
    plot_snp_distribution(df, OUTPUT_DIR)
    
    print("--- Script Finished Successfully ---")
