import re
import pandas as pd
import matplotlib.pyplot as plt

# Input and output
log_file = "verlet_debug.log"
csv_file = "verlet_debug.csv"

# Regex to capture the debug lines
pattern = re.compile(
    r"\[DEBUG_VERLET\] Energy change in step: ([\-\d\.eE]+)% \(TE_pre=([\-\d\.eE]+), TE_post=([\-\d\.eE]+)\)"
)

data = []
with open(log_file, "r") as f:
    for line in f:
        match = pattern.search(line)
        if match:
            drift = float(match.group(1))
            TE_pre = float(match.group(2))
            TE_post = float(match.group(3))
            data.append([drift, TE_pre, TE_post])

# Convert to DataFrame
df = pd.DataFrame(data, columns=["drift_percent", "TE_pre", "TE_post"])
df["step"] = range(len(df))

# Save to CSV
df.to_csv(csv_file, index=False)
print(f"Saved {len(df)} entries to {csv_file}")

# Plot drift
plt.figure(figsize=(10,6))
plt.plot(df["step"], df["drift_percent"], marker="o", linestyle="-", label="Energy Drift %")
plt.axhline(0, color="black", linestyle="--")
plt.xlabel("Step")
plt.ylabel("Energy Drift (%)")
plt.title("Energy Drift per Verlet Step")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("verlet_drift_plot.png")
plt.show()
