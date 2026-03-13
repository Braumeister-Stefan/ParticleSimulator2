import pandas as pd
import matplotlib.pyplot as plt

# Load your simulation data
df = pd.read_csv('Inputs/rendered_scenarios/TestCase2.csv')

print("Data overview:")
print(f"Total steps: {df['step_id'].max() + 1}")
print(f"Total particles: {df[' particle_id'].nunique()}")

# Check if positions change
particle_0 = df[df[' particle_id'] == 0]
print(f"\nParticle 0 position range:")
print(f"X: {particle_0['x'].min():.6f} to {particle_0['x'].max():.6f}")
print(f"Y: {particle_0['y'].min():.6f} to {particle_0['y'].max():.6f}")

# Check multiple particles
for pid in df[' particle_id'].unique()[:3]:  # Check first 3 particles
    particle_data = df[df[' particle_id'] == pid]
    x_range = particle_data['x'].max() - particle_data['x'].min()
    y_range = particle_data['y'].max() - particle_data['y'].min()
    print(f"Particle {pid}: X movement={x_range:.6f}, Y movement={y_range:.6f}")

# Quick plot of trajectories
plt.figure(figsize=(10, 8))
for pid in df[' particle_id'].unique():
    particle_data = df[df[' particle_id'] == pid]
    plt.plot(particle_data['x'], particle_data['y'], label=f'Particle {pid}', linewidth=2)
    plt.scatter(particle_data['x'].iloc[0], particle_data['y'].iloc[0], 
                marker='o', s=100, color='red')  # Start position

plt.xlabel('X Position')
plt.ylabel('Y Position')
plt.title('Particle Trajectories - TestCase2')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.savefig('trajectory_check.png', dpi=150, bbox_inches='tight')
plt.show()

print("\nIf trajectories are just points, particles aren't moving!")