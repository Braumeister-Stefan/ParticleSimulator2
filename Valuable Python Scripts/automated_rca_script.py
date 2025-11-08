import os
import pandas as pd
import numpy as np

# Paths
input_folder = "Inputs/rendered_scenarios"
test_cases = {
    "MVP Test (2 particles)": "MVP_Test_metrics.csv",
    "Validation Test (10 particles)": "TestCase2_metrics.csv",
    "Inelastic Test": "Inelastic_Test_metrics.csv",
    "LargeDT Test": "LargeDT_Test_metrics.csv"
}

def analyze_te_drift(df, label):
    initial_TE = df["TE"].iloc[0]
    rel_error = (df["TE"] - initial_TE) / abs(initial_TE) * 100
    max_drift = rel_error.abs().max()
    
    # Identify step where TE jump is largest
    te_diff = df["TE"].diff().fillna(0)
    max_jump_step = df["step_id"].iloc[te_diff.abs().idxmax()]
    max_jump_value = te_diff.abs().max()
    
    # Basic automated suggestion
    suggestion = ""
    if max_drift > 50:  # huge drift
        suggestion = "Check multi-particle force calculations or integration stability"
    elif max_drift > 1:
        suggestion = "Minor drift: possibly timestep effect (Δt)"
    else:
        suggestion = "Stable: no major energy conservation issues"
    
    return {
        "Test Case": label,
        "Max TE Drift (%)": max_drift,
        "Step of Max TE Jump": max_jump_step,
        "Max TE Jump": max_jump_value,
        "Suggested Focus": suggestion
    }

# Collect results
results = []

for label, filename in test_cases.items():
    file_path = os.path.join(input_folder, filename)
    if not os.path.exists(file_path):
        print(f"{filename} not found, skipping {label}.")
        continue
    df = pd.read_csv(file_path)
    if "TE" not in df.columns or "step_id" not in df.columns:
        print(f"{filename} missing required columns, skipping {label}.")
        continue
    analysis = analyze_te_drift(df, label)
    results.append(analysis)

# Convert to DataFrame for easy view
results_df = pd.DataFrame(results)
print("\n=== Automated Energy Drift Analysis ===")
print(results_df)

# Optional: save summary CSV
results_df.to_csv("Plots/Automated_RCA_Summary.csv", index=False)
print("\nSummary saved as 'Plots/Automated_RCA_Summary.csv'")
