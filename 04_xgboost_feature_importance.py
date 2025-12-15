#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script Name: 04_xgboost_feature_importance.py
Description:
    This script performs feature importance analysis using XGBoost.
    It ranks phenotypic traits based on their contribution to the target trait (Lodging Score).
    Features are grouped (Morphological, Mechanical, etc.) for better visualization.

    Includes synthetic data generation for privacy protection and reproducibility.

Author:      Liang Xiaotian
Email:       494382219@qq.com
Date:        2025/12/06
License:     MIT
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.preprocessing import StandardScaler

# Try to import xgboost, handle error if missing
try:
    import xgboost as xgb
except ImportError:
    print("[Error] 'xgboost' library not found. Please install it using: pip install xgboost")
    sys.exit(1)

# ==========================================
# 1. Configuration & Setup
# ==========================================

# Set random seed for reproducibility
np.random.seed(2025)

# Create output directory (Consistent with R scripts)
output_dir = "results"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Set font for plots (Robust handling)
try:
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
except Exception:
    pass # Fallback to default
plt.rcParams['axes.unicode_minus'] = False

# ==========================================
# 2. Data Preparation (Synthetic Data)
# ==========================================

def generate_synthetic_data(n_samples=300):
    """
    Generates synthetic phenotypic data for demonstration purposes.
    Simulates the structure of oat lodging traits.
    """
    print("[Data] Generating synthetic data...")
    
    # Define feature groups based on the manuscript
    # Group names aligned with biological meaning
    groups = {
        "Morphological": ["PH", "PL", "CGH", "SCIL", "SCID"],
        "Basal_1st":     ["FIL", "FID", "FIWT", "FIBR", "FIPS"],
        "Basal_2nd":     ["SIL", "SID", "SIWT", "SIBR", "SIPS"],
        "Basal_3rd":     ["TIL", "TID", "TIWT", "TIBR", "TIPS"]
    }
    
    # Flatten feature list
    features = [f for g in groups.values() for f in g]
    
    # Generate random feature data
    X = pd.DataFrame(np.random.normal(loc=10, scale=2, size=(n_samples, len(features))), columns=features)
    
    # Generate target variable (Lodging Score, LS)
    # Simulate: LS is negatively correlated with stem strength (TIBR) and positively with Height (PH)
    noise = np.random.normal(0, 0.5, n_samples)
    y = 0.6 * X['PH'] - 0.4 * X['TIBR'] - 0.3 * X['SIBR'] + 0.1 * X['TID'] + noise
    
    # Normalize y to 0-3 range (Lodging Score)
    y = ((y - y.min()) / (y.max() - y.min())) * 3
    
    # Combine into one dataframe for consistency
    data = X.copy()
    data['LS'] = y
    data.insert(0, 'SampleID', [f'Sample_{i+1}' for i in range(n_samples)])
    
    return data, groups

# Load data
data, feature_groups = generate_synthetic_data()

# ==========================================
# 3. Data Preprocessing
# ==========================================

print("[Processing] Preprocessing data...")

# Separate features (X) and target (y)
target_col = 'LS'
y = data[target_col]
X = data.drop(columns=['SampleID', target_col])

# Standardization (Important for comparing importance correctly)
scaler = StandardScaler()
X_scaled = pd.DataFrame(scaler.fit_transform(X), columns=X.columns)

# ==========================================
# 4. Model Training (XGBoost)
# ==========================================

print(f"[Analysis] Training XGBoost model to predict {target_col}...")

# Hyperparameters from the manuscript (Section 2.3.2)
xgb_model = xgb.XGBRegressor(
    objective="reg:squarederror",
    max_depth=3,
    learning_rate=0.1,
    n_estimators=100,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=2025,
    n_jobs=-1
)

xgb_model.fit(X_scaled, y)

# ==========================================
# 5. Feature Importance Extraction
# ==========================================

print("[Analysis] Extracting and grouping feature importance...")

# Extract raw importance scores
importance_df = pd.DataFrame({
    "Feature": X_scaled.columns,
    "Score": xgb_model.feature_importances_
})

# Assign groups to features
def assign_group(feature_name):
    for group, members in feature_groups.items():
        if feature_name in members:
            return group
    return "Other"

importance_df["Group"] = importance_df["Feature"].apply(assign_group)

# Define group order for plotting
group_order = ["Morphological", "Basal_1st", "Basal_2nd", "Basal_3rd", "Other"]
importance_df["Group"] = pd.Categorical(importance_df["Group"], categories=group_order, ordered=True)

# Sort by Group and then by Score
importance_sorted = importance_df.sort_values(by=["Group", "Score"], ascending=[True, False]).reset_index(drop=True)

# Save results
csv_path = os.path.join(output_dir, "04_XGBoost_Feature_Importance.csv")
importance_sorted.to_csv(csv_path, index=False)
print(f"    - CSV saved to: {csv_path}")

# ==========================================
# 6. Visualization
# ==========================================

print("[Plotting] Generating feature importance plot...")

# Define color palette (Publication-ready)
group_colors = {
    "Morphological": "#1f77b4", # Blue
    "Basal_1st":     "#ff7f0e", # Orange
    "Basal_2nd":     "#2ca02c", # Green
    "Basal_3rd":     "#d62728", # Red
    "Other":         "#9467bd"  # Purple
}

plt.figure(figsize=(10, 8))

# Draw Horizontal Bars
for idx, row in importance_sorted.iterrows():
    color = group_colors.get(row["Group"], "#333333")
    plt.barh(y=idx, width=row["Score"], color=color, height=0.7, align='center')
    
    # Add value labels
    plt.text(row["Score"] + 0.002, idx, f'{row["Score"]:.3f}', va='center', fontsize=8)

# Formatting
plt.yticks(range(len(importance_sorted)), importance_sorted["Feature"])
plt.xlabel("Importance Score (Gain)", fontsize=12, fontweight='bold')
plt.ylabel("Traits", fontsize=12, fontweight='bold')
plt.title(f"Feature Importance for Lodging Resistance (XGBoost)", fontsize=14, fontweight='bold')

# Invert Y-axis (Top importance at top)
plt.gca().invert_yaxis()

# Create Legend
patches = [Patch(color=color, label=group) for group, color in group_colors.items() if group in importance_sorted["Group"].unique()]
plt.legend(handles=patches, title="Trait Group", loc="lower right")

plt.tight_layout()

# Save Plot
plot_path = os.path.join(output_dir, "04_XGBoost_Feature_Importance.png")
plt.savefig(plot_path, dpi=300)
print(f"    - Plot saved to: {plot_path}")

print("\n[Success] File 04 Analysis Complete.")
