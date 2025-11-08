import pandas as pd
import matplotlib.pyplot as plt

# File path for Validation Test (10 particles)
file = "Inputs/rendered_scenarios/TestCase2_metrics.csv"

# Load metrics
df = pd.read_csv(file)

# Rename columns if needed (your file uses step_id not step)
if "step_id" in df.columns:
    df.rename(columns={"step_id": "step"}, inplace=True)

# Zoom window
start, end = 700, 900
df_window = df[(df["step"] >= start) & (df["step"] <= end)]

# Plot KE, PE, TE
plt.figure(figsize=(12,6))
plt.plot(df_window["step"], df_window["KE"], label="Kinetic Energy (KE)")
plt.plot(df_window["step"], df_window["PE"], label="Potential Energy (PE)")
plt.plot(df_window["step"], df_window["TE"], label="Total Energy (TE)")
plt.xlabel("Step")
plt.ylabel("Energy")
plt.title(f"Zoomed Energies (Steps {start}-{end})")
plt.legend()
plt.grid(True)
plt.show()

# Plot Momentum (x, y)
plt.figure(figsize=(12,6))
plt.plot(df_window["step"], df_window["Px"], label="Momentum X")
plt.plot(df_window["step"], df_window["Py"], label="Momentum Y")
plt.xlabel("Step")
plt.ylabel("Momentum")
plt.title(f"Zoomed Momentum (Steps {start}-{end})")
plt.legend()
plt.grid(True)
plt.show()
