import os
import pandas as pd
import numpy as np
from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.models import LinearAxis, Range1d, FuncTickFormatter

# ========== USER CONFIGURABLE PARAMETERS ==========
G = 6.67430e-11
csv_folder = r"C:\Users\smdw1\OneDrive\Bureaublad\Development\Projects\cpp\ParticleSimulator2\Inputs\rendered_scenarios"

THRESHOLD_VALUE = 100
ADDITIONAL_STEPS_PERCENT = 0.5

redraw_plot = True
zoom_focus_step = 6400
zoom_window = 25

ZOOM_SELECTIONS = ["Error", "y=0", "TE"]
# ========== END OF USER CONFIGURABLE PARAMETERS ==========

def find_stable_cutoff(energies, threshold):
    length = len(energies)
    for i in range(length):
        if energies[i:].max() - energies[i:].min() < threshold:
            return i
    return length - 1

def compute_energies(df):
    def kinetic_energy(subdf):
        if 'vz' in subdf.columns:
            return 0.5 * subdf['m'] * (subdf['vx']**2 + subdf['vy']**2 + subdf['vz']**2)
        else:
            return 0.5 * subdf['m'] * (subdf['vx']**2 + subdf['vy']**2)
    def potential_energy(group):
        x = group['x'].values
        y = group['y'].values
        z = group['z'].values if 'z' in group.columns else np.zeros(len(group))
        m = group['m'].values
        n = len(m)
        if n < 2:
            return pd.Series(0.0, index=group.index)
        pe_sum = 0.0
        for i in range(n):
            for j in range(i+1, n):
                dx = x[i] - x[j]
                dy = y[i] - y[j]
                dz = z[i] - z[j]
                r = np.sqrt(dx*dx + dy*dy + dz*dz)
                pe_sum += -G*m[i]*m[j]/r
        return pd.Series(pe_sum / n, index=group.index)
    df['kinetic_energy'] = kinetic_energy(df)
    if 'step_id' not in df.columns:
        df['step_id'] = 0
    df['potential_energy'] = df.groupby('step_id', group_keys=False).apply(potential_energy)
    df['heating_energy'] = df['temp'] if 'temp' in df.columns else 0.0
    df['total_energy'] = df['kinetic_energy'] + df['potential_energy'] + df['heating_energy']
    return df

def compute_error(actual_energy, initial_energy):
    return actual_energy - initial_energy

def create_TE_and_error_plot(x, TE, marginal_error, title, TE0, err_range,
                             te_min, te_max, width=500, height=400,
                             cpp_series=None, te_label="Total Energy (py)",
                             error_label="Marginal Error", y0_label="y=0",
                             plot_TE=True):
    p = figure(title=title, x_axis_label='Step Index', width=width, height=height)
    p.y_range = Range1d(start=te_min, end=te_max)
    p.yaxis.axis_label = "Energy (J)"

    if plot_TE:
        if te_label is not None:
            p.line(x, TE, color="black", line_width=2, legend_label=te_label)
        else:
            p.line(x, TE, color="black", line_width=2)

    if cpp_series is not None:
        p.line(x, cpp_series, color="gray", line_width=2, line_dash="dashed",
               legend_label="Total Energy (cpp)")

    transformed_error = marginal_error + TE0
    p.extra_y_ranges = {"error": Range1d(start=TE0 - err_range, end=TE0 + err_range)}
    formatter_code = f"return (tick - {TE0}).toFixed(2);"
    right_axis = LinearAxis(
        y_range_name="error",
        axis_label="Marginal Error (J)",
        formatter=FuncTickFormatter(code=formatter_code)
    )
    p.add_layout(right_axis, 'right')

    if error_label is not None:
        p.line(x, transformed_error, color="red", line_dash="dotted", line_width=2,
               legend_label=error_label, y_range_name="error")
    else:
        p.line(x, transformed_error, color="red", line_dash="dotted", line_width=2,
               y_range_name="error")

    if y0_label is not None:
        p.line(x, [TE0]*len(x), color="gray", line_width=1,
               legend_label=y0_label, y_range_name="error")
    else:
        p.line(x, [TE0]*len(x), color="gray", line_width=1, y_range_name="error")

    p.legend.location = "top_right"
    return p

def plot_results(base_name, final_energy_df, stable_step_id):
    x = final_energy_df.index.values
    TE = final_energy_df['total_energy'].values
    KE = final_energy_df['kinetic_energy'].values
    PE = final_energy_df['potential_energy'].values
    HE = final_energy_df['heating_energy'].values

    ERR_cumulative = final_energy_df['error'].values
    marginal_error = np.diff(ERR_cumulative, prepend=ERR_cumulative[0])
    marginal_error[0] = 0

    TE0 = TE[0]
    err_range = np.max(np.abs(marginal_error))
    te_min = np.min(TE)
    te_max = np.max(TE)

    cpp_series = final_energy_df["total_energy_cpp"].values if "total_energy_cpp" in final_energy_df.columns else None

    # (1) Full-range TE + Marginal Error
    p1 = create_TE_and_error_plot(
        x, TE, marginal_error,
        title=f"{base_name} - Total Energy & Marginal Error",
        TE0=TE0, err_range=err_range, te_min=te_min, te_max=te_max,
        cpp_series=cpp_series,
        te_label="Total Energy (py)",
        error_label="Marginal Error",
        y0_label="y=0",
        plot_TE=True
    )

    if "margin_TE_error" in final_energy_df.columns:
        margin_cpp = final_energy_df["margin_TE_error"].values
        p1.line(x, margin_cpp, color="gray", line_dash="dotted",
                legend_label="Marginal Error (cpp)", y_range_name="error")
        p1.scatter(x, margin_cpp, marker="asterisk", size=6,
                   color="gray", legend_label="Marginal Error (cpp)", y_range_name="error")

    # (2) TE Zoom
    zoom_start = max(0, zoom_focus_step - zoom_window)
    zoom_end = zoom_focus_step + zoom_window
    zoom_mask = (x >= zoom_start) & (x <= zoom_end)

    if not np.any(zoom_mask):
        x_zoom = x
        TE_zoom = TE
        marginal_error_zoom = marginal_error
    else:
        x_zoom = x[zoom_mask]
        TE_zoom = TE[zoom_mask]
        marginal_error_zoom = marginal_error[zoom_mask]

    zoom_te_label = "Total Energy (py)" if "TE" in ZOOM_SELECTIONS else None
    zoom_error_label = "Marginal Error" if "Error" in ZOOM_SELECTIONS else None
    zoom_y0_label = "y=0" if "y=0" in ZOOM_SELECTIONS else None
    plot_TE_zoom = True if "TE" in ZOOM_SELECTIONS else False

    p2 = create_TE_and_error_plot(
        x_zoom, TE_zoom, marginal_error_zoom,
        title=f"Zoom around step={zoom_focus_step} ± {zoom_window}",
        TE0=TE0, err_range=err_range, te_min=te_min, te_max=te_max,
        cpp_series=None,
        te_label=zoom_te_label,
        error_label=zoom_error_label,
        y0_label=zoom_y0_label,
        plot_TE=plot_TE_zoom
    )

    # (3) Kinetic, Potential, & Heating
    p3 = figure(title="Kinetic, Potential, & Heating",
                x_axis_label="Step Index", y_axis_label="Energy (J)",
                width=500, height=400)
    p3.line(x, KE, color="blue", line_width=2, legend_label="Kinetic")
    p3.line(x, PE, color="red", line_width=2, legend_label="Potential")
    p3.line(x, HE, color="green", line_width=2, legend_label="Heating")
    p3.legend.location = "top_right"

    # (4) Overlap Marginal Errors (bar for transl & corrected, line for overlap_iters)
    p4 = figure(title="Overlap Marginal Errors",
                x_axis_label="Step Index", y_axis_label="Margin (J)",
                width=500, height=400)
    x_list = list(x)

    # We want separate y-axis for "overlap_iters_in_step":
    p4.extra_y_ranges = {"iters": Range1d(start=0, end=10)}  # We'll adjust dynamically
    p4.add_layout(LinearAxis(y_range_name="iters", axis_label="Overlap Iters"), 'right')

    # (A) If columns exist, bar-plot them side by side
    bar_width = 0.4

    def bar_if_exists(colname, color, offset):
        if colname in final_energy_df.columns:
            vals = final_energy_df[colname].values
            x_offset = [xx + offset for xx in x_list]
            p4.vbar(x=x_offset, top=vals, width=bar_width, color=color,
                    legend_label=colname)

    # margin_TE_error_overlap_ij_transl & margin_TE_error_overlap_ij_corrected
    bar_if_exists("margin_TE_error_overlap_ij_transl", "blue", 0.0)
    bar_if_exists("margin_TE_error_overlap_ij_corrected", "red", 0.4)

    # (B) If "overlap_iters_in_step" exists, line-plot it w/ dotted style on right axis
    if "overlap_iters_in_step" in final_energy_df.columns:
        iters_vals = final_energy_df["overlap_iters_in_step"].values
        p4.extra_y_ranges["iters"].end = max(1, max(iters_vals) * 1.2)
        p4.line(x_list, iters_vals, y_range_name="iters",
                line_dash="dotted", color="gray", line_width=2,
                legend_label="overlap_iters_in_step")

    p4.legend.location = "top_right"

    # (5) Zoom for Overlap Marginal Errors
    p5 = figure(title=f"Overlap Marginal Errors (Zoom) {zoom_focus_step} ± {zoom_window}",
                x_axis_label="Step Index", y_axis_label="Margin (J)",
                width=500, height=400)

    # Similar approach for the zoom
    x_list_zoom = x_list if not np.any(zoom_mask) else [xx for (xx, m) in zip(x_list, zoom_mask) if m]

    p5.extra_y_ranges = {"iters": Range1d(start=0, end=10)}
    p5.add_layout(LinearAxis(y_range_name="iters", axis_label="Overlap Iters"), 'right')

    def bar_zoom_if_exists(colname, color, offset):
        if colname in final_energy_df.columns:
            full_vals = final_energy_df[colname].values
            if not np.any(zoom_mask):
                zoom_vals = full_vals
            else:
                zoom_vals = [val for (val, m) in zip(full_vals, zoom_mask) if m]
            x_offset_zoom = [xx + offset for xx in x_list_zoom]
            p5.vbar(x=x_offset_zoom, top=zoom_vals, width=bar_width, color=color,
                    legend_label=colname)

    # Bars
    bar_zoom_if_exists("margin_TE_error_overlap_ij_transl", "blue", 0.0)
    bar_zoom_if_exists("margin_TE_error_overlap_ij_corrected", "red", 0.4)

    # Overlap iters line
    if "overlap_iters_in_step" in final_energy_df.columns:
        full_iters = final_energy_df["overlap_iters_in_step"].values
        if not np.any(zoom_mask):
            iters_zoom = full_iters
        else:
            iters_zoom = [val for (val, m) in zip(full_iters, zoom_mask) if m]
        p5.extra_y_ranges["iters"].end = max(1, max(iters_zoom) * 1.2)
        p5.line(x_list_zoom, iters_zoom, y_range_name="iters",
                line_dash="dotted", color="gray", line_width=2,
                legend_label="overlap_iters_in_step")

    p5.legend.location = "top_right"

    from bokeh.layouts import gridplot
    grid = gridplot([
        [p1, p2],
        [p3, p4],
        [p5, None]
    ])

    output_path = os.path.join(csv_folder, f"{base_name}_Enhanced_Energy_Visualization.html")
    output_file(output_path, title=f"{base_name} Enhanced Energy Visualization")
    save(grid)
    print(f"Generated Bokeh figure for '{base_name}' saved to {output_path}")

def main():
    processed_files = 0
    skipped_files = 0
    for filename in os.listdir(csv_folder):
        if filename.lower().endswith(".csv") and "_metrics" not in filename.lower():
            csv_path = os.path.join(csv_folder, filename)
            base_name = os.path.splitext(filename)[0]
            output_file_path = os.path.join(csv_folder, f"{base_name}_Enhanced_Energy_Visualization.html")

            if os.path.isfile(output_file_path) and not redraw_plot:
                print(f"Skipping '{filename}' - plot already exists.")
                skipped_files += 1
                continue

            print(f"Processing file: {filename}")
            df = pd.read_csv(csv_path)
            df = compute_energies(df)

            energy_df = df.groupby('step_id', as_index=True)[
                ['kinetic_energy', 'potential_energy', 'heating_energy', 'total_energy']
            ].sum().sort_index()

            total_energy_series = energy_df['total_energy'].values
            stable_position = find_stable_cutoff(total_energy_series, THRESHOLD_VALUE)
            total_steps = len(energy_df.index)

            additional_steps = int(ADDITIONAL_STEPS_PERCENT * stable_position)
            cutoff_position = stable_position + additional_steps
            if cutoff_position >= total_steps:
                cutoff_position = total_steps - 1

            step_ids = energy_df.index.tolist()
            stable_step_id = step_ids[cutoff_position]
            df = df[df['step_id'] <= stable_step_id].copy()

            final_energy_df = df.groupby('step_id', as_index=True)[
                ['kinetic_energy', 'potential_energy', 'heating_energy', 'total_energy']
            ].sum().sort_index()

            metrics_filename = base_name + "_metrics.csv"
            metrics_path = os.path.join(csv_folder, metrics_filename)
            if os.path.isfile(metrics_path):
                try:
                    metrics_df = pd.read_csv(metrics_path)
                    metrics_df = metrics_df.set_index("step_id")

                    if "TE" in metrics_df.columns:
                        final_energy_df = final_energy_df.join(
                            metrics_df["TE"].rename("total_energy_cpp"), how="left")

                    # Additional possible columns
                    for col in metrics_df.columns:
                        if col not in final_energy_df.columns:
                            final_energy_df[col] = metrics_df[col]

                    print(f"Associated metrics file found: {metrics_filename}")
                except Exception as e:
                    print(f"Error processing metrics file {metrics_filename}: {e}")

            # Replace infinities with NaN, then fill any NaN with 0.0
            final_energy_df.replace([np.inf, -np.inf], np.nan, inplace=True)
            final_energy_df.fillna(0.0, inplace=True)

            initial_total_energy = final_energy_df['total_energy'].iloc[0]
            final_energy_df['error'] = compute_error(
                final_energy_df['total_energy'], initial_total_energy)

            plot_results(base_name, final_energy_df, stable_step_id)
            processed_files += 1

    print("===== Run Summary =====")
    print(f"Processed files: {processed_files}")
    print(f"Skipped files:   {skipped_files}")

if __name__ == "__main__":
    if False:
        main()
