import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Gravitational constant (m^3 kg^-1 s^-2)
G = 6.67430e-11

# Folder containing CSV files
csv_folder = r"C:\Users\smdw1\OneDrive\Bureaublad\Development\Projects\cpp\ParticleSimulator2\Inputs\rendered_scenarios"

# User-defined parameters
THRESHOLD_VALUE = 10   # Range threshold to consider system stable
redraw_plot = True     # Regenerate plots even if they exist

processed_files = 0
skipped_files = 0

def find_stable_cutoff(energies, threshold):
    """
    Identify the earliest step index at which the remaining
    portion of energies (from that index to the end) stays
    within a 'threshold' range (max-min < threshold).
    If no such index is found, returns the last index.
    """
    length = len(energies)
    for i in range(length):
        slice_range = energies[i:].max() - energies[i:].min()
        if slice_range < threshold:
            return i
    return length - 1

# Loop through each CSV in the folder
for filename in os.listdir(csv_folder):
    if filename.lower().endswith(".csv"):
        csv_path = os.path.join(csv_folder, filename)
        base_name = os.path.splitext(filename)[0]
        energy_plot_path = os.path.join(csv_folder, f"{base_name}_Enhanced_Energy_Visualization.jpg")

        # Skip if output already exists (and not forcing redraw)
        if os.path.isfile(energy_plot_path) and not redraw_plot:
            print(f"Skipping '{filename}' - plot already exists.")
            skipped_files += 1
            continue

        print(f"Processing file: {filename}")
        df = pd.read_csv(csv_path)

        # Remove columns z, vz if present
        columns_to_drop = ['z', 'vz']
        df.drop(columns=[col for col in columns_to_drop if col in df.columns],
                errors='ignore', inplace=True)

        # Check for required columns
        required_cols = ['x', 'y', 'vx', 'vy', 'm']
        if not all(col in df.columns for col in required_cols):
            print(f"Skipping '{filename}' - missing required columns.")
            skipped_files += 1
            continue

        # 1. Kinetic Energy
        df['kinetic_energy'] = 0.5 * df['m'] * (df['vx']**2 + df['vy']**2)

        # 2. Potential Energy
        def calculate_potential(group):
            x = group['x'].values
            y = group['y'].values
            m = group['m'].values
            n = len(m)
            if n < 2:
                return pd.Series(0.0, index=group.index)

            dx = x[:, None] - x
            dy = y[:, None] - y
            r = np.sqrt(dx**2 + dy**2)
            np.fill_diagonal(r, np.inf)  # Avoid division by zero

            m_matrix = m[:, None] * m
            pe_terms = G * m_matrix / r
            return pd.Series(-pe_terms.sum(axis=1), index=group.index)

        df['potential_energy'] = df.groupby('step_id', group_keys=False).apply(calculate_potential)

        # 3. Heating Energy (from 'temp' if available, else 0)
        if 'temp' in df.columns:
            df['heating_energy'] = df['temp']
        else:
            df['heating_energy'] = 0.0

        # 4. Total Energy = KE + PE + HE
        df['total_energy'] = df['kinetic_energy'] + df['potential_energy'] + df['heating_energy']

        # 5. Group by step & find cutoff
        energy_df = df.groupby('step_id', as_index=True)[
            ['kinetic_energy', 'potential_energy', 'heating_energy', 'total_energy']
        ].mean().sort_index()

        step_ids = energy_df.index.to_list()
        total_energy_series = energy_df['total_energy'].values

        stable_position = find_stable_cutoff(total_energy_series, THRESHOLD_VALUE)
        total_steps = len(step_ids)
        additional_steps = int(0.1 * stable_position)
        cutoff_position = stable_position + additional_steps
        if cutoff_position >= total_steps:
            cutoff_position = total_steps - 1

        stable_step_id = step_ids[cutoff_position]
        df = df[df['step_id'] <= stable_step_id].copy()

        # Recompute energies after truncation
        final_energy_df = df.groupby('step_id', as_index=True)[
            ['kinetic_energy', 'potential_energy', 'heating_energy', 'total_energy']
        ].mean().sort_index()

        # Prepare data for plotting
        step_index = final_energy_df.index
        KE = final_energy_df['kinetic_energy']
        PE = final_energy_df['potential_energy']
        HE = final_energy_df['heating_energy']
        TE = final_energy_df['total_energy']

        # For the kinetic & potential subplot, show P - K:
        p_minus_k = PE + KE

        # For the kinetic & heating subplot, show sum (K + H):
        k_plus_h = KE + HE

        # -----------------------------------------------------------------
        # Create a figure with FOUR subplots:
        #
        #  1) Top-most:       Total Energy alone
        #  2) Next:           Kinetic, Potential, Heating
        #  3) Next:           Kinetic vs. Potential + difference (P-K)
        #  4) Bottom-most:    Kinetic vs. Heating + sum (K+H)
        # -----------------------------------------------------------------
        fig, axes = plt.subplots(4, 1, figsize=(12, 18), sharex=True)
        ax_top      = axes[0]
        ax_bottom   = axes[1]
        ax_kp       = axes[2]
        ax_kh       = axes[3]

        # (1) Plot Total Energy
        ax_top.plot(step_index, TE, color='black', label='Total Energy')
        ax_top.set_title(f"{base_name} - Total Energy")
        ax_top.set_ylabel('Energy (J)')
        ax_top.legend(loc='upper right')
        ax_top.grid(True, alpha=0.3)

        # (2) Plot Kinetic, Potential, Heating
        ax_bottom.plot(step_index, KE, color='blue', label='Kinetic')
        ax_bottom.plot(step_index, PE, color='red', label='Potential')
        ax_bottom.plot(step_index, HE, color='green', label='Heating')
        ax_bottom.set_title("Kinetic, Potential, & Heating")
        ax_bottom.set_ylabel('Energy (J)')
        ax_bottom.legend(loc='upper right')
        ax_bottom.grid(True, alpha=0.3)

        # (3) Kinetic vs. Potential, plus (P - K)
        ax_kp.plot(step_index, KE, color='blue', label='Kinetic (K)')
        ax_kp.plot(step_index, PE, color='red',  label='Potential (P)')
        ax_kp.plot(step_index, p_minus_k, color='gray', linestyle='--',
                   label='Sum (P + K)')
        ax_kp.set_title("Kinetic vs. Potential + Sum (P + K)")
        ax_kp.set_ylabel('Energy (J)')
        ax_kp.legend(loc='upper right')
        ax_kp.grid(True, alpha=0.3)

        # (4) Kinetic vs. Heating, plus (K + H)
        ax_kh.plot(step_index, KE, color='blue', label='Kinetic (K)')
        ax_kh.plot(step_index, HE, color='green', label='Heating (H)')
        ax_kh.plot(step_index, k_plus_h, color='gray', linestyle='--',
                   label='Sum (K + H)')
        ax_kh.set_title("Kinetic vs. Heating + Sum (K + H)")
        ax_kh.set_xlabel('Step ID')
        ax_kh.set_ylabel('Energy (J)')
        ax_kh.legend(loc='upper right')
        ax_kh.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(energy_plot_path, dpi=300, bbox_inches='tight')
        plt.close()

        processed_files += 1
        print(f"Generated enhanced 4-subplot energy plot for '{filename}'.\n")

# Final summary
print("===== Run Summary =====")
print(f"Processed files: {processed_files}")
print(f"Skipped files:   {skipped_files}")
