import pandas as pd
import matplotlib.pyplot as plt

# Load metrics
df = pd.read_csv("Inputs/rendered_scenarios/TestCase2_metrics.csv")

# Focus on critical window
window = df[(df["step_id"] >= 700) & (df["step_id"] <= 830)]

# Print some stats
print(window[["step_id", "KE", "PE", "TE"]])

# Plot
plt.figure(figsize=(10,6))
plt.plot(window["step_id"], window["KE"], label="Kinetic Energy")
plt.plot(window["step_id"], window["PE"], label="Potential Energy")
plt.plot(window["step_id"], window["TE"], label="Total Energy")
plt.xlabel("Step")
plt.ylabel("Energy")
plt.title("Energy Dynamics (steps 700–830)")
plt.legend()
plt.grid(True)
plt.show()
