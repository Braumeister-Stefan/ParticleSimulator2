import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- Step 1: Load / Create Drift DataFrame ---
# Example: simulate data for 100 particles and 100 steps (replace this with your real log data)
num_particles = 100
num_steps = 100

# Example: random drift data for demonstration (replace with your actual data)
np.random.seed(0)
drift_data = np.random.normal(loc=0, scale=1e-6, size=(num_steps, num_particles))
drift_df = pd.DataFrame(drift_data, columns=[f'particle_{i}' for i in range(num_particles)])
drift_df['step'] = range(num_steps)

# --- Step 2: Compute per-particle statistics ---
particle_stats = drift_df.drop(columns='step').agg(['mean', 'std', 'min', 'max']).T
print("Per-Particle Drift Statistics:")
print(particle_stats.head())

# --- Step 3: Detect anomalies ---
# Define threshold for unusual drift (adjust based on expected range)
threshold = -1e-5
anomalies = (drift_df.drop(columns='step') < threshold)
anomaly_counts = anomalies.sum(axis=0)
print("\nParticles with anomalies (count of steps exceeding threshold):")
print(anomaly_counts[anomaly_counts > 0])

# --- Step 4: Heatmap of per-particle drift over steps ---
plt.figure(figsize=(15, 8))
sns.heatmap(drift_df.drop(columns='step').T, cmap='coolwarm', cbar_kws={'label': 'Drift'})
plt.xlabel("Step")
plt.ylabel("Particle")
plt.title("Heatmap of Per-Particle Drift Over Steps")
plt.show()

# --- Step 5: Highlight high-drift particles ---
high_drift_particles = anomaly_counts[anomaly_counts > 0].index.tolist()
print("\nHigh-drift particles detected:", high_drift_particles)
