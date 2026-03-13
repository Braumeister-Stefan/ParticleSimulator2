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

# Initialize summary
summary = []

# Loop through test cases
for label, filename in test_cases.items():
    file_path = os.path.join(input_folder, filename)
    if not os.path.exists(file_path):
        print(f"Warning: {filename} not found. Skipping {label}.")
        continue

    df = pd.read_csv(file_path)

    # Check required columns
    required_cols = ["step_id", "KE", "PE", "TE"]
    if not all(col in df.columns for col in required_cols):
        print(f"Skipping {filename}: required columns missing")
        continue

    # Compute relative TE error
    initial_TE = df["TE"].iloc[0]
    rel_error = (df["TE"] - initial_TE)/abs(initial_TE)*100
    max_drift = rel_error.abs().max()
    summary.append((label, max_drift))

    # 1️⃣ Energies vs Step
    plt.figure(figsize=(10,6))
    plt.plot(df["step_id"], df["KE"], label="KE", color="blue")
    plt.plot(df["step_id"], df["PE"], label="PE", color="green")
    plt.plot(df["step_id"], df["TE"], label="TE", color="red")
    plt.xlabel("Step")
    plt.ylabel("Energy")
    plt.title(f"Energies vs Step: {label}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{label.replace(' ', '_')}_energies.png"))
    plt.close()

    # 2️⃣ Relative TE Error
    plt.figure(figsize=(10,6))
    plt.plot(df["step_id"], rel_error, label="Relative TE Error (%)", color="purple")
    plt.axhline(0, color="black", linestyle="--", linewidth=1)
    plt.axhline(1, color="red", linestyle="--", linewidth=1, label="±1% Threshold")
    plt.axhline(-1, color="red", linestyle="--", linewidth=1)
    plt.xlabel("Step")
    plt.ylabel("Relative TE Error (%)")
    plt.title(f"Total Energy Drift: {label}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{label.replace(' ', '_')}_TE_error.png"))
    plt.close()

    # 3️⃣ Optional: Momentum plots if columns exist
    if "momentum_x" in df.columns and "momentum_y" in df.columns:
        plt.figure(figsize=(10,6))
        plt.plot(df["step_id"], df["momentum_x"], label="Momentum X", color="orange")
        plt.plot(df["step_id"], df["momentum_y"], label="Momentum Y", color="cyan")
        plt.xlabel("Step")
        plt.ylabel("Momentum")
        plt.title(f"Momentum vs Step: {label}")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, f"{label.replace(' ', '_')}_momentum.png"))
        plt.close()

# Print summary table
print("\nMaximum TE drift per test case:")
for case, drift in summary:
    print(f"{case}: {drift:.3f}%")

print("\nAll plots have been saved in the 'Plots' folder.")
