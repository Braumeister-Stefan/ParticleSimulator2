import pandas as pd
import matplotlib.pyplot as plt
import re

# Read log file
with open("verlet_steps.log") as f:
    lines = f.readlines()

pattern = re.compile(r"\[DEBUG_VERLET_STEP\] step=(\d+), drift=([-0-9.eE]+)%,")

data = []
for line in lines:
    match = pattern.search(line)
    if match:
        step = int(match.group(1))
        drift = float(match.group(2))
        data.append((step, drift))

df = pd.DataFrame(data, columns=["step", "drift_percent"])
df.to_csv("verlet_drift_debug.csv", index=False)
print(df.head())

# Plot
plt.figure(figsize=(8,5))
plt.plot(df["step"], df["drift_percent"], marker="o", label="Drift %")
plt.axhline(0, color="black", linestyle="--")
plt.xlabel("Step")
plt.ylabel("Relative TE Drift (%)")
plt.title("Verlet Step Energy Drift (First 100 Steps)")
plt.legend()
plt.grid(True)
plt.savefig("verlet_drift_debug.png", dpi=150)
plt.show()
