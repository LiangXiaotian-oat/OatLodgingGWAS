#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SNP Genotype Visualization Script
---------------------------
Copyright (c) 2025 Liang Xiaotian <494382219@qq.com>
--------------------------------------------------
Description: 
    This script reads SNP genotype data from a local file (Excel/CSV) 
    and generates a grouped bar chart visualization.
    It automatically handles dependencies and saves outputs to a specific folder.

Usage:
    Ensure your data file (e.g., 'snp_data.xlsx') is in the same directory.
    Run the script: python snp_visualization.py

Author: Liang Xiaotian
Email: 494382219@qq.com
Date: 2025/12/09
License: MIT
"""

import os
import sys
import subprocess

# ==========================================
# 1. Dependency Management & Auto-Install
# ==========================================
def install_package(package):
    """Install a package using pip if import fails."""
    try:
        print(f"Installing missing package: {package}...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    except subprocess.CalledProcessError:
        print(f"Error: Failed to install {package}. Please install it manually.")
        sys.exit(1)

# List of required packages
required_packages = {
    'pandas': 'pandas',
    'matplotlib': 'matplotlib',
    'seaborn': 'seaborn',
    'openpyxl': 'openpyxl' # Required for reading Excel files
}

# Check and import packages dynamically
for lib_name, pip_name in required_packages.items():
    try:
        __import__(lib_name)
    except ImportError:
        install_package(pip_name)

# Safe imports after check
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# ==========================================
# 2. Configuration (User Inputs)
# ==========================================

# File Settings
# INPUT_FILE: The name of your local Excel or CSV file.
INPUT_FILE = "single_snp_分组柱状图.xlsx" 

# Output Settings
OUTPUT_DIR = "results"
FIG_SIZE = (12, 7)
DPI_SETTING = 300

# Font Settings
# Try to use 'Arial' or standard sans-serif
FONT_NAME = 'Arial' 

# Color Settings (Academic Style)
COLOR_PALETTE = {
    "A": "#D6221E",  # Red
    "T": "#3778AD",  # Blue
    "C": "#4EA74A",  # Green
    "G": "#EF7C1B",  # Orange/Yellow
    "N": "gray"      # Missing/Null
}

# ==========================================
# 3. Data Loading Logic
# ==========================================
def load_dataset(file_path):
    """
    Loads data from Excel or CSV file.
    """
    if not os.path.exists(file_path):
        print(f"[Error] File not found: {file_path}")
        print("Please ensure the data file is in the current directory.")
        sys.exit(1)

    try:
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        elif file_path.endswith(('.xls', '.xlsx')):
            df = pd.read_excel(file_path)
        else:
            print("[Error] Unsupported file format. Please use .csv or .xlsx")
            sys.exit(1)
        
        print(f"Successfully loaded data from: {file_path}")
        return df
    except Exception as e:
        print(f"[Error] Failed to read file: {e}")
        sys.exit(1)

# ==========================================
# 4. Plotting Function
# ==========================================
def plot_snp_distribution(df, output_folder):
    """
    Plots the grouped bar chart and saves it to the output folder.
    """
    # Create folder if not exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output directory: {output_folder}")

    # Set up the plot style
    sns.set_style("ticks")
    plt.figure(figsize=FIG_SIZE)
    
    # Font configuration
    plt.rcParams['font.family'] = FONT_NAME
    plt.rcParams['font.size'] = 12

    # Draw the barplot
    try:
        ax = sns.barplot(
            data=df,
            x="SNP",
            y="Count",
            hue="Genotype",
            palette=COLOR_PALETTE,
            edgecolor="black",
            linewidth=1.2
        )
    except ValueError as e:
        print(f"[Error] Plotting failed. Check column names (SNP, Count, Genotype). Details: {e}")
        sys.exit(1)

    # Labels and Title
    plt.title("Single SNP Genotype Counts", fontsize=16, fontweight='bold', pad=20)
    plt.xlabel("SNP Position", fontsize=14, fontweight='bold')
    plt.ylabel("Sample Count", fontsize=14, fontweight='bold')

    # Adjust ticks
    plt.xticks(rotation=45, ha="right")

    # Legend optimization (Outside the plot area)
    plt.legend(title="Genotype", title_fontsize='12', 
               bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)

    # Add value labels on top of bars
    for container in ax.containers:
        ax.bar_label(container, padding=2, fontsize=10)

    plt.tight_layout()

    # Save logic
    # 1. Save as PNG (Raster)
    png_path = os.path.join(output_folder, "SNP_Genotype_Analysis.png")
    plt.savefig(png_path, dpi=DPI_SETTING, bbox_inches='tight')
    
    # 2. Save as PDF (Vector)
    pdf_path = os.path.join(output_folder, "SNP_Genotype_Analysis.pdf")
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')

    print(f"Success! Figures saved to:\n - {png_path}\n - {pdf_path}")
    
    # Display plot
    plt.show()

# ==========================================
# 5. Main Execution
# ==========================================
if __name__ == "__main__":
    print("--- Starting SNP Visualization Script ---")
    
    # 1. Load Data
    df = load_dataset(INPUT_FILE)
    
    # 2. Preview Data
    print("Data Preview (First 5 rows):")
    print(df.head())
    
    # 3. Generate Visualization
    plot_snp_distribution(df, OUTPUT_DIR)
    
    print("--- Script Finished Successfully ---")
