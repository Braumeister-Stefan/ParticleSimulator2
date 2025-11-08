import os
import pandas as pd
import matplotlib.pyplot as plt

# Paths
input_folder = "Inputs/rendered_scenarios"
output_folder = "Plots"
os.makedirs(output_folder, exist_ok=True)

# Test case metrics CSVs
test_cases = {
    "MVP Test (2 particles)": "MVP_Test_metrics.csv",
    "Validation Test (10 particles)": "TestCase2_metrics.csv",
    "Inelastic Test": "Inelastic_Test_metrics.csv",
    "LargeDT Test": "LargeDT_Test_metrics.csv"
}

plt.figure(figsize=(12,7))

for label, filename in test_cases.items():
    file_path = os.path.join(input_folder, filename)
    if not os.path.exists(file_path):
        continue

    df = pd.read_csv(file_path)
    if "TE" not in df.columns or "step_id" not in df.columns:
        continue

    initial_TE = df["TE"].iloc[0]
    rel_error = (df["TE"] - initial_TE)/abs(initial_TE)*100

    plt.plot(df["step_id"], rel_error, label=f"{label} (max {rel_error.abs().max():.2f}%)")

    # Highlight step with max TE jump for Validation Test
    if label == "Validation Test (10 particles)":
        te_diff = df["TE"].diff().fillna(0)
        max_jump_idx = te_diff.abs().idxmax()
        max_jump_step = df["step_id"].iloc[max_jump_idx]
        max_jump_value = rel_error.iloc[max_jump_idx]
        plt.scatter(max_jump_step, max_jump_value, color="red", s=100, zorder=5,
                    label=f"TE Blow-up Step {max_jump_step}")

# Reference lines
plt.axhline(0, color="black", linestyle="--", linewidth=1)
plt.axhline(1, color="red", linestyle="--", linewidth=1, label="±1% Threshold")
plt.axhline(-1, color="red", linestyle="--", linewidth=1)

plt.xlabel("Step")
plt.ylabel("Relative TE Error (%)")
plt.title("Total Energy Drift Comparison with Highlighted Blow-up")
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save combined figure
combined_plot_path = os.path.join(output_folder, "Combined_TE_drift_with_blowup.png")
plt.savefig(combined_plot_path)
plt.show()

print(f"Combined TE drift plot with blow-up saved as: {combined_plot_path}")
