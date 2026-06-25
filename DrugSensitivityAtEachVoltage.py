import os
import numpy as np
import pandas as pd

# Size of each "cell" in the TIFF heatmap (pixels per drug and per voltage)
BLOCK_SIZE = 100

# For heatmap TIFF
try:
    import tifffile as tiff
    HAVE_TIFFFILE = True
except ImportError:
    HAVE_TIFFFILE = False
    print("WARNING: 'tifffile' not found. Heatmap TIFF will NOT be created.")
    print("Install it with:  pip install tifffile")

# For labeled PNG with colorbar
try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False
    print("WARNING: 'matplotlib' not found. Labeled PNG heatmap will NOT be created.")
    print("Install it with:  pip install matplotlib")


# -------------------------------------------------------
# Utility functions
# -------------------------------------------------------

def safe_sem(values):
    """
    Compute SEM (standard error of the mean) ignoring NaNs.
    Returns NaN if fewer than 2 non-NaN values.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[~np.isnan(arr)]
    n = arr.size
    if n <= 1:
        return np.nan
    return float(np.nanstd(arr, ddof=1) / np.sqrt(n))


def sem_of_block_percent(mean_ctrl, sem_ctrl, mean_drug, sem_drug):
    """
    Approximate SEM for block% using error propagation.

    block% = 100 * (1 - mean_drug / mean_ctrl)
    """
    if (
        mean_ctrl is None
        or mean_drug is None
        or np.isnan(mean_ctrl)
        or np.isnan(mean_drug)
        or mean_ctrl == 0
    ):
        return np.nan

    var_mc = sem_ctrl ** 2 if not np.isnan(sem_ctrl) else 0.0
    var_md = sem_drug ** 2 if not np.isnan(sem_drug) else 0.0

    # f = 100 * (1 - md / mc)
    # ∂f/∂mc = 100 * (md / mc^2)
    # ∂f/∂md = -100 * (1 / mc)
    df_dmc = 100.0 * (mean_drug / (mean_ctrl ** 2))
    df_dmd = -100.0 * (1.0 / mean_ctrl)

    var_block = (df_dmc ** 2) * var_mc + (df_dmd ** 2) * var_md
    return float(np.sqrt(var_block))


def split_control_drug_columns(df, genotype_prefix, voltage_col="Voltage (mV)"):
    """
    Split columns into control vs drug for a given genotype (e.g. 'wt' or 'dKO').

    Heuristic as before: classify based on stripped tokens of the column name.
    """
    prefix_lower = genotype_prefix.lower()
    control_cols = []
    drug_cols = []

    for col in df.columns:
        if col == voltage_col:
            continue
        s = str(col).strip()
        if not s.lower().startswith(prefix_lower):
            continue

        tokens = s.split()
        if len(tokens) == 1:
            control_cols.append(col)
        else:
            if tokens[1].startswith("."):
                control_cols.append(col)
            else:
                drug_cols.append(col)

    return control_cols, drug_cols


def stats_for_group(df, row_idx, cols):
    """
    For a single voltage row (row_idx) and a list of columns (cols),
    return (n, mean_abs_current, sem_abs_current).
    Uses absolute value of current.
    """
    if not cols:
        return 0, np.nan, np.nan

    try:
        vals = df.loc[row_idx, cols].to_numpy(dtype=float)
    except Exception:
        vals = []
        for c in cols:
            try:
                vals.append(float(df.loc[row_idx, c]))
            except Exception:
                vals.append(np.nan)
        vals = np.array(vals, dtype=float)

    vals = np.abs(vals)
    vals = vals[~np.isnan(vals)]
    n = vals.size
    if n == 0:
        return 0, np.nan, np.nan

    mean_val = float(np.nanmean(vals))
    sem_val = safe_sem(vals)
    return n, mean_val, sem_val


def create_fire_ice_lut():
    """
    Fire/Ice LUT WITHOUT white at the hot end.
    - Cold = cyan → blue
    - Neutral = black (0)
    - Hot = dark red → red → bright orange/yellow
    """
    # Control points in [0,1]
    xp = np.array([0.0, 0.25, 0.5, 0.75, 0.9, 1.0], dtype=float)

    # Colors as RGB 0-255:
    #   cyan, blue, black, dark red, red, bright yellow/orange
    cp = np.array([
        [0,   255, 255],   # cyan
        [0,   0,   255],   # blue
        [0,   0,   0],     # black
        [139, 0,   0],     # dark red
        [255, 0,   0],     # red
        [255, 200, 0],     # bright yellow/orange
    ], dtype=float)

    lut = np.zeros((256, 3), dtype=np.uint8)
    x = np.linspace(0, 1, 256)
    for i in range(3):
        lut[:, i] = np.interp(x, xp, cp[:, i]).astype(np.uint8)
    return lut


def create_fire_ice_cmap():
    """
    Matplotlib Fire/Ice colormap — no white.
    Negative = icy (cyan→blue), 0 = black, Positive = red→orange→yellow.
    """
    colors = [
        (0.0,  (0/255,   1.0,     1.0)),      # cyan
        (0.25, (0/255,   0/255,   1.0)),      # blue
        (0.50, (0.0,     0.0,     0.0)),      # black
        (0.75, (139/255, 0.0,     0.0)),      # dark red
        (0.90, (1.0,     0.0,     0.0)),      # red
        (1.00, (1.0,     200/255, 0.0)),      # bright yellow/orange
    ]
    return LinearSegmentedColormap.from_list("fire_ice_nowhite", colors)


# -------------------------------------------------------
# File processing: per-voltage drug sensitivity
# -------------------------------------------------------

def process_file(path, v_min, v_max, voltage_col="Voltage (mV)"):
    """
    Process one CSV file to compute drug block at EACH voltage.

    Returns:
      per_voltage_rows: list of dicts (for CSV)
      voltage_to_ratio: { voltage_mV : sensitivity_ratio_WT_over_dKO }
    """
    fname = os.path.basename(path)
    df = pd.read_csv(path)

    if voltage_col not in df.columns:
        raise ValueError(f"{fname}: Voltage column '{voltage_col}' not found.")

    wt_ctrl_cols, wt_drug_cols = split_control_drug_columns(df, "wt", voltage_col)
    dko_ctrl_cols, dko_drug_cols = split_control_drug_columns(df, "dKO", voltage_col)

    print(f"\nFile: {fname}")
    print(f"  WT control cols: {wt_ctrl_cols}")
    print(f"  WT drug cols:    {wt_drug_cols}")
    print(f"  dKO control cols: {dko_ctrl_cols}")
    print(f"  dKO drug cols:    {dko_drug_cols}")

    voltages = df[voltage_col].to_numpy(dtype=float)

    # Mask for voltage range
    mask = np.ones_like(voltages, dtype=bool)
    if v_min is not None:
        mask &= (voltages >= v_min)
    if v_max is not None:
        mask &= (voltages <= v_max)

    per_voltage_rows = []
    voltage_to_ratio = {}

    for idx in range(len(df)):
        if not mask[idx]:
            continue

        v = float(voltages[idx])

        wt_n_ctrl, wt_mean_ctrl, wt_sem_ctrl = stats_for_group(df, idx, wt_ctrl_cols)
        wt_n_drug, wt_mean_drug, wt_sem_drug = stats_for_group(df, idx, wt_drug_cols)

        dko_n_ctrl, dko_mean_ctrl, dko_sem_ctrl = stats_for_group(df, idx, dko_ctrl_cols)
        dko_n_drug, dko_mean_drug, dko_sem_drug = stats_for_group(df, idx, dko_drug_cols)

        def block_percent(mean_ctrl, mean_drug):
            if (
                mean_ctrl is None
                or np.isnan(mean_ctrl)
                or mean_ctrl == 0
                or np.isnan(mean_drug)
            ):
                return np.nan
            return float(100.0 * (1.0 - (mean_drug / mean_ctrl)))

        wt_block = block_percent(wt_mean_ctrl, wt_mean_drug)
        dko_block = block_percent(dko_mean_ctrl, dko_mean_drug)

        wt_block_sem = sem_of_block_percent(
            wt_mean_ctrl, wt_sem_ctrl, wt_mean_drug, wt_sem_drug
        )
        dko_block_sem = sem_of_block_percent(
            dko_mean_ctrl, dko_sem_ctrl, dko_mean_drug, dko_sem_drug
        )

        # ---- NEW: sensitivity ratio uses the MAGNITUDE of block (abs), so it is never negative ----
        # This avoids cases where tiny apparent facilitation in dKO (negative block)
        # would yield a negative WT/dKO ratio, which is not biologically what you mean by "sensitivity".
        if (
            wt_block is None
            or dko_block is None
            or np.isnan(wt_block)
            or np.isnan(dko_block)
            or dko_block == 0
        ):
            sens_ratio = np.nan
        else:
            sens_ratio = float(abs(wt_block) / abs(dko_block))

        voltage_to_ratio[v] = sens_ratio

        row = {
            "filename": fname,
            "Voltage_mV": v,

            "wt_n_control": wt_n_ctrl,
            "wt_n_drug": wt_n_drug,
            "dko_n_control": dko_n_ctrl,
            "dko_n_drug": dko_n_drug,

            "wt_mean_control_Current": wt_mean_ctrl,
            "wt_sem_control_Current": wt_sem_ctrl,
            "wt_mean_drug_Current": wt_mean_drug,
            "wt_sem_drug_Current": wt_sem_drug,

            "wt_block_percent_mean": wt_block,
            "wt_block_percent_sem": wt_block_sem,

            "dko_mean_control_Current": dko_mean_ctrl,
            "dko_sem_control_Current": dko_sem_ctrl,
            "dko_mean_drug_Current": dko_mean_drug,
            "dko_sem_drug_Current": dko_sem_drug,

            "dko_block_percent_mean": dko_block,
            "dko_block_percent_sem": dko_block_sem,

            "sensitivity_ratio_WT_over_dKO": sens_ratio,
        }

        per_voltage_rows.append(row)

    return per_voltage_rows, voltage_to_ratio


# -------------------------------------------------------
# Main
# -------------------------------------------------------

def main():
    # 1. Directory with CSV files
    dir_path = input("Enter the path to the directory containing the .csv files: ").strip()
    if not os.path.isdir(dir_path):
        print(f"ERROR: '{dir_path}' is not a valid directory.")
        return

    # 2. Voltage range
    print("\nSpecify the voltage range (in mV) over which to analyze drug sensitivity.")
    print("Press Enter to use the FULL voltage range in each file.\n")

    v_min = None
    v_max = None

    v_min_str = input("Enter minimum voltage (mV) (or press Enter for full range): ").strip()
    if v_min_str != "":
        try:
            v_min = float(v_min_str)
        except ValueError:
            print(f"WARNING: Could not parse '{v_min_str}' as float; using no minimum (full range).")
            v_min = None

    v_max_str = input("Enter maximum voltage (mV) (or press Enter for full range): ").strip()
    if v_max_str != "":
        try:
            v_max = float(v_max_str)
        except ValueError:
            print(f"WARNING: Could not parse '{v_max_str}' as float; using no maximum (full range).")
            v_max = None

    print(
        f"\nUsing voltage range: "
        f"{v_min if v_min is not None else 'min(data)'} to "
        f"{v_max if v_max is not None else 'max(data)'} mV\n"
    )

    # 3. Collect CSV files
    csv_files = [f for f in os.listdir(dir_path) if f.lower().endswith(".csv")]
    csv_files.sort()

    if not csv_files:
        print("No .csv files found in the specified directory.")
        return

    all_rows = []
    heatmap_info = []  # list of (filename, {voltage: sensitivity_ratio})

    # 4. Process each file
    for fname in csv_files:
        full_path = os.path.join(dir_path, fname)
        try:
            rows, v2ratio = process_file(
                full_path, v_min, v_max, voltage_col="Voltage (mV)"
            )
            all_rows.extend(rows)
            heatmap_info.append((fname, v2ratio))
        except Exception as e:
            print(f"\nERROR while processing '{fname}': {e}")
            continue

    # 5. Write per-voltage summary CSV
    if all_rows:
        summary_df = pd.DataFrame(all_rows)
        out_csv = os.path.join(dir_path, "Drug sensitivity at each voltage.csv")
        summary_df.to_csv(out_csv, index=False)
        print("\nPer-voltage drug sensitivity written to:")
        print(f"  {out_csv}")
    else:
        print("No per-voltage rows generated (check for errors above).")

    # 6. Build heatmap data
    if not heatmap_info:
        print("\nNo heatmap data collected.")
        return

    all_voltages = set()
    for _, v2ratio in heatmap_info:
        all_voltages.update(v2ratio.keys())

    if not all_voltages:
        print("\nNo voltages available for heatmap.")
        return

    voltages_sorted = sorted(all_voltages)
    n_volt = len(voltages_sorted)
    n_drugs = len(heatmap_info)

    heatmap_matrix = np.full((n_drugs, n_volt), np.nan, dtype=np.float32)
    drug_labels = []

    for row_idx, (fname, v2ratio) in enumerate(heatmap_info):
        drug_labels.append(fname)
        for col_idx, v in enumerate(voltages_sorted):
            val = v2ratio.get(v, np.nan)
            heatmap_matrix[row_idx, col_idx] = np.float32(val if val is not None else np.nan)

    # ----- Determine symmetric scale around 0 for visualization -----
    valid = np.isfinite(heatmap_matrix)
    if np.any(valid):
        data_min = float(np.nanmin(heatmap_matrix[valid]))
        data_max = float(np.nanmax(heatmap_matrix[valid]))
        max_abs = max(abs(data_min), abs(data_max))
        if max_abs == 0:
            max_abs = 1.0
    else:
        data_min, data_max, max_abs = 0.0, 0.0, 1.0

    vmin = -max_abs
    vmax = max_abs

    # -------- 8-bit scaling for TIFF: -max_abs -> 0, 0 -> 128, +max_abs -> 255 --------
    scaled = np.zeros_like(heatmap_matrix, dtype=np.uint8)
    if np.any(valid):
        norm = (heatmap_matrix[valid] + max_abs) / (2 * max_abs)  # -max_abs->0, 0->0.5, +max_abs->1
        norm = np.clip(norm, 0.0, 1.0)
        scaled_vals = (255.0 * norm).astype(np.uint8)
        scaled[valid] = scaled_vals
    # NaN entries remain 0 (cold side)

    # -------- Enlarge each cell to BLOCK_SIZE x BLOCK_SIZE and save TIFF --------
    if HAVE_TIFFFILE:
        if BLOCK_SIZE > 1:
            scaled_big = np.repeat(np.repeat(scaled, BLOCK_SIZE, axis=0), BLOCK_SIZE, axis=1)
        else:
            scaled_big = scaled

        heatmap_path = os.path.join(dir_path, "Drug_sensitivity_heatmap_8bit.tif")
        try:
            tiff.imwrite(heatmap_path, scaled_big)
            print("\n8-bit HEATMAP written to:")
            print(f"  {heatmap_path}")
            print("Open in ImageJ.")
            print("  - Pixel intensity encodes WT/dKO block sensitivity ratio magnitude")
            print("  - 0     → cold extreme (~ -max ratio)")
            print("  - 128   → 0 (black, no bias)")
            print("  - 255   → hot extreme (~ +max ratio)")
            print("  - X-axis (columns) = Voltage (mV)")
            print("  - Y-axis (rows)    = File/Drug")
            print(f"  - Each cell is {BLOCK_SIZE}x{BLOCK_SIZE} pixels")
        except Exception as e:
            print(f"\nERROR writing heatmap TIFF: {e}")
    else:
        print("\nSkipping TIFF heatmap because tifffile is not installed.")

    # -------- Write Fire/Ice LUT for ImageJ --------
    try:
        lut = create_fire_ice_lut()
        lut_path = os.path.join(dir_path, "Drug_sensitivity_FireIce.lut")
        lut.astype(np.uint8).tofile(lut_path)
        print("\nFire/Ice LUT written to:")
        print(f"  {lut_path}")
        print("In ImageJ/Fiji: open the TIFF, then Image > Color > Load LUT... and select this file.")
    except Exception as e:
        print(f"\nERROR writing Fire/Ice LUT file: {e}")

    # Axis label CSVs
    try:
        volts_csv = os.path.join(dir_path, "Drug_sensitivity_heatmap_voltages.csv")
        drugs_csv = os.path.join(dir_path, "Drug_sensitivity_heatmap_drugs.csv")

        pd.DataFrame({"Voltage_mV": voltages_sorted}).to_csv(volts_csv, index=False)
        pd.DataFrame({"row_index": np.arange(n_drugs),
                      "filename": drug_labels}).to_csv(drugs_csv, index=False)

        print("\nHeatmap axis label CSVs written to:")
        print(f"  {volts_csv}")
        print(f"  {drugs_csv}")
    except Exception as e:
        print(f"\nERROR writing heatmap label CSVs: {e}")

    # -------- Labeled PNG with axes + colorbar (Fire/Ice colormap) --------
    if HAVE_MPL:
        try:
            fire_ice_cmap = create_fire_ice_cmap()

            fig, ax = plt.subplots(figsize=(max(6, n_volt * 0.3), max(4, n_drugs * 0.3)))

            if np.any(valid):
                im = ax.imshow(
                    heatmap_matrix,
                    aspect="auto",
                    interpolation="nearest",
                    cmap=fire_ice_cmap,
                    vmin=vmin,
                    vmax=vmax,
                )
            else:
                im = ax.imshow(
                    heatmap_matrix,
                    aspect="auto",
                    interpolation="nearest",
                    cmap=fire_ice_cmap,
                )

            # X ticks = voltages (subsample if many)
            if n_volt <= 20:
                xticks = np.arange(n_volt)
            else:
                xticks = np.linspace(0, n_volt - 1, num=10, dtype=int)
            ax.set_xticks(xticks)
            ax.set_xticklabels([f"{voltages_sorted[i]:.0f}" for i in xticks], rotation=90)

            # Y ticks = each drug/file
            ax.set_yticks(np.arange(n_drugs))
            ax.set_yticklabels(drug_labels)

            ax.set_xlabel("Voltage (mV)")
            ax.set_ylabel("Drug / File")
            ax.set_title("WT/dKO Block Sensitivity Ratio (|WT block| / |dKO block|)")

            cbar = fig.colorbar(im, ax=ax)
            cbar.set_label("WT/dKO block sensitivity ratio (magnitude)", rotation=270, labelpad=15)

            if np.any(valid):
                cbar_ticks = [vmin, 0.0, vmax]
                cbar.set_ticks(cbar_ticks)
                cbar.set_ticklabels([f"{vmin:.1f}", "0", f"{vmax:.1f}"])

            fig.tight_layout()

            png_path = os.path.join(dir_path, "Drug_sensitivity_heatmap_labeled.png")
            fig.savefig(png_path, dpi=300)
            plt.close(fig)

            print("\nLabeled PNG heatmap written to:")
            print(f"  {png_path}")
            print("This image has:")
            print("  - X-axis labeled with Voltage (mV)")
            print("  - Y-axis labeled with Drug/File names")
            print("  - Fire/Ice colorbar centered at 0 (black)")
            print("    • cold (icy) low ratio")
            print("    • hot (fiery) high ratio")
        except Exception as e:
            print(f"\nERROR creating labeled PNG heatmap: {e}")
    else:
        print("\nSkipping labeled PNG because matplotlib is not installed.")


if __name__ == "__main__":
    main()
