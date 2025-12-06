#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
XGBoost Feature Importance Analysis
-----------------------------------
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
    This script demonstrates how to use XGBoost to analyze the importance of 
    phenotypic traits in predicting lodging resistance (or other target traits).
    It includes data simulation, model training, feature grouping, and 
    visualization of feature importance.

Usage:
    python xgboost_feature_analysis.py

Author: Liang Xiaotian
Email: 494382219@qq.com
Date: 2025/12/06
"""

import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

# ==========================================
# 1. Configuration & Setup
# ==========================================

# Set random seed for reproducibility
np.random.seed(42)

# Create output directory
output_dir = "XGBoost_Analysis_Output"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Set font for plots (Arial is standard for scientific publications)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial'] 
# If you need Chinese characters, uncomment the line below and ensure SimHei is installed
# plt.rcParams['font.sans-serif'] = ['SimHei'] 
plt.rcParams['axes.unicode_minus'] = False

# ==========================================
# 2. Data Preparation (Synthetic Data)
# ==========================================

def generate_synthetic_data(n_samples=300):
    """
    Generates synthetic phenotypic data for demonstration purposes.
    In real usage, replace this with pd.read_csv('your_data.csv').
    """
    print("[*] Generating synthetic data...")
    
    # Define feature groups
    groups = {
        "Group1": ["PH", "SL", "SFW", "GCH", "RD", "RL"],
        "Group2": ["FIL", "FID", "FIBR", "FIPS", "FIWT", "FIFW", "FIDW"],
        "Group3": ["SIL", "SID", "SIBR", "SIPS", "SIWT", "SIFW", "SIDW"],
        "Group4": ["TIL", "TID", "TIBR", "TIPS", "TIWT", "TIFW", "TIDW"]
    }
    
    # Flatten feature list
    features = [f for g in groups.values() for f in g]
    
    # Generate random data
    data = pd.DataFrame(np.random.randn(n_samples, len(features)), columns=features)
    
    # Generate target variable (LD) with some dependency on features
    # Let's say PH (Plant Height) and SIBR (Stem Strength) are important
    noise = np.random.normal(0, 0.5, n_samples)
    data['LD'] = 0.5 * data['PH'] - 0.3 * data['SIBR'] + 0.1 * data['TID'] + noise
    
    # Add ID column
    data.insert(0, 'ID', [f'Sample_{i}' for i in range(n_samples)])
    
    return data, groups

# Load data (Using synthetic data here)
# To use your own data: 
# data = pd.read_csv('wj+cz初筛.csv')
# feature_groups = { ... define your groups ... }
data, feature_groups = generate_synthetic_data()

# ==========================================
# 3. Data Preprocessing
# ==========================================

print("[*] Preprocessing data...")

# Check for missing values
if data.isnull().sum().sum() > 0:
    print("Warning: Missing values detected. Filling with mean.")
    data = data.fillna(data.mean(numeric_only=True))

# Separate features (X) and target (y)
y = data['LD']
X = data.drop(columns=['ID', 'LD'])

# Standardization
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
X_scaled = pd.DataFrame(X_scaled, columns=X.columns)

# ==========================================
# 4. Model Training (XGBoost)
# ==========================================

print("[*] Training XGBoost model...")

xgb_model = xgb.XGBRegressor(
    objective="reg:squarederror",
    max_depth=3,
    learning_rate=0.1,
    n_estimators=100,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=42,
    n_jobs=-1
)

xgb_model.fit(X_scaled, y)

# ==========================================
# 5. Feature Importance Analysis
# ==========================================

print("[*] Extracting feature importance...")

# Extract raw importance scores
feature_importance = pd.DataFrame({
    "Feature": X_scaled.columns,
    "Score": xgb_model.feature_importances_
})

# Assign groups to features
def get_feature_group(feature):
    for group_name, feats in feature_groups.items():
        if feature in feats:
            return group_name
    return "Other"

feature_importance["Group"] = feature_importance["Feature"].apply(get_feature_group)

# Sort by Group Order then by Score Descending
group_order = list(feature_groups.keys()) + ["Other"]
feature_importance["Group"] = pd.Categorical(
    feature_importance["Group"],
    categories=group_order,
    ordered=True
)

feature_importance_sorted = feature_importance.sort_values(
    by=["Group", "Score"],
    ascending=[True, False]
).reset_index(drop=True)

# Export results table
csv_path = os.path.join(output_dir, "Feature_Importance_Sorted.csv")
feature_importance_sorted.to_csv(csv_path, index=False)
print(f"    - Results saved to: {csv_path}")

# ==========================================
# 6. Visualization
# ==========================================

print("[*] Generating visualization...")

# Define colors for groups
# Using a colorblind-friendly palette
group_colors = {
    "Group1": "#1f77b4",  # Blue
    "Group2": "#ff7f0e",  # Orange
    "Group3": "#2ca02c",  # Green
    "Group4": "#d62728",  # Red
    "Other":  "#9467bd"   # Purple
}

plt.figure(figsize=(12, 10))

# Plot bars
for idx, row in feature_importance_sorted.iterrows():
    feature = row["Feature"]
    score = row["Score"]
    group = row["Group"]
    color = group_colors.get(group, group_colors["Other"])
    
    # Horizontal bar
    plt.barh(y=idx, width=score, color=color, height=0.7)
    
    # Add text label (5 decimal places)
    plt.text(
        x=score + (max(feature_importance_sorted['Score']) * 0.01), # Dynamic offset
        y=idx,
        s=f"{score:.5f}",
        ha="left",
        va="center",
        fontsize=9
    )

# Configure axes
plt.yticks(
    ticks=range(len(feature_importance_sorted)),
    labels=feature_importance_sorted["Feature"],
    fontsize=10
)
plt.xlabel("Importance Score", fontsize=12)
plt.ylabel("Trait / Feature", fontsize=12)
plt.title("XGBoost Feature Importance (Grouped)", fontsize=14, pad=20)

# Invert Y axis to have the first group at the top
plt.gca().invert_yaxis()

# Add Legend
legend_elements = [
    Patch(facecolor=color, edgecolor='none', label=group)
    for group, color in group_colors.items()
    if group in feature_importance_sorted["Group"].unique()
]
plt.legend(handles=legend_elements, title="Feature Group", loc="lower right", fontsize=10)

plt.tight_layout()

# Save plot
plot_path = os.path.join(output_dir, "XGBoost_Feature_Importance.png")
plt.savefig(plot_path, dpi=300, bbox_inches="tight")
plt.close()

print(f"    - Plot saved to: {plot_path}")
print("\n[Success] Analysis complete.")
