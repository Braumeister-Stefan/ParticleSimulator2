import pandas as pd
import re

# --- Step 1: Define log file ---
log_file = "verlet_steps.log"
output_csv = "verlet_parsed.csv"

# --- Step 2: Define regex pattern ---
# Adjust the regex based on how your log records appear.
# Example pattern assumes lines like:
# [DEBUG_VERLET_STEP] step=123, particle=5, drift=1.2e-06, force=0.005, work=0.0001, distance=0.02
pattern = re.compile(
    r"\[DEBUG_VERLET_STEP\]\s+step=(\d+),\s+particle=(\d+),\s+drift=([-0-9.eE]+),\s+force=([-0-9.eE]+),\s+work=([-0-9.eE]+),\s+distance=([-0-9.eE]+)"
)

# --- Step 3: Parse the log file ---
data = []
with open(log_file, "r") as f:
    for line in f:
        match = pattern.search(line)
        if match:
            step = int(match.group(1))
            particle = int(match.group(2))
            drift = float(match.group(3))
            force = float(match.group(4))
            work = float(match.group(5))
            distance = float(match.group(6))
            data.append([step, particle, drift, force, work, distance])

# --- Step 4: Convert to DataFrame ---
columns = ["step", "particle", "drift", "force", "work", "distance"]
df = pd.DataFrame(data, columns=columns)

# --- Step 5: Save to CSV ---
df.to_csv(output_csv, index=False)
print(f"Parsed log saved to {output_csv}, total records: {len(df)}")
