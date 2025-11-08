import pandas as pd
import matplotlib.pyplot as plt
import re

# --- Step 1: Read the log file ---
log_file = "verlet_steps.log"  # adjust path if needed
with open(log_file, "r") as f:
    lines = f.readlines()

# --- Step 2: Define regex pattern ---
# Matches lines like: [DEBUG_VERLET_STEP] step=123, drift=0.456%
pattern = re.compile(r"\[DEBUG_VERLET_STEP\] step=(\d+), drift=([-0-9.eE]+)%")

# --- Step 3: Extract data ---
data = []
for line in lines:
    match = pattern.search(line)
    if match:
        step = int(match.group(1))
        drift = float(match.group(2))
        data.append([step, drift])

# --- Step 4: Convert to DataFrame ---
df = pd.DataFrame(data, columns=["Step", "Drift"])

# --- Step 5: Plot the drift vs. steps ---
plt.figure(figsize=(10, 6))
plt.plot(df["Step"], df["Drift"], marker='o', linestyle='-', color='blue', label='Drift')
plt.xlabel("Step")
plt.ylabel("Drift (%)")
plt.title("Verlet Steps vs. Drift")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
