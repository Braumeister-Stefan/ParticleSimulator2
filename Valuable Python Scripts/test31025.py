import pandas as pd
import re
import matplotlib.pyplot as plt

# =============================
# Step 1: Load / parse log
# =============================
log_file = "verlet_steps.log"
output_csv = "verlet_energy_log.csv"

# Regex for current log format: step, drift, TE_pre, TE_post
pattern = re.compile(
    r"\[DEBUG_VERLET_STEP\]\s+step=(\d+),\s+drift=([-0-9.eE]+)%,\s+TE_pre=([-0-9.eE]+),\s+TE_post=([-0-9.eE]+)"
)

data = []
with open(log_file, "r") as f:
    for line in f:
        match = pattern.search(line)
        if match:
            step = int(match.group(1))
            drift = float(match.group(2))
            TE_pre = float(match.group(3))
            TE_post = float(match.group(4))
            data.append([step, drift, TE_pre, TE_post])

# Create DataFrame
df = pd.DataFrame(data, columns=["step", "drift", "TE_pre", "TE_post"])
df["delta_TE"] = df["TE_post"] - df["TE_pre"]

# Save CSV for later use
df.to_csv(output_csv, index=False)
print(f"Parsed log saved to {output_csv}, total records: {len(df)}")

# =============================
# Step 2: Detect unstable steps
# =============================
threshold = df["delta_TE"].abs().mean() + 3 * df["delta_TE"].abs().std()
unstable_steps = df[df["delta_TE"].abs() > threshold]

print(f"\nUnstable steps detected: {len(unstable_steps)}")
print(unstable_steps[["step", "drift", "TE_pre", "TE_post", "delta_TE"]])

# =============================
# Step 3: Plot Total Energy vs Time
# =============================
plt.figure(figsize=(12,6))
plt.plot(df["step"], df["TE_post"], label="TE_post", color="blue")
plt.plot(df["step"], df["TE_pre"], label="TE_pre", color="orange", alpha=0.7)
plt.xlabel("Simulation Step")
plt.ylabel("Total Energy (TE)")
plt.title("Total Energy vs Simulation Steps")
plt.legend()
plt.grid(True)

# Highlight unstable steps on the plot
plt.scatter(unstable_steps["step"], unstable_steps["TE_post"], color="red", label="Unstable Steps")
plt.legend()
plt.savefig("Plots/TE_vs_time.png", dpi=300)
plt.show()

# =============================
# Step 4: Summary Table for Report
# =============================
summary_file = "Plots/unstable_steps_summary.csv"
unstable_steps.to_csv(summary_file, index=False)
print(f"Summary of unstable steps saved to {summary_file}")
