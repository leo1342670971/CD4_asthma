import argparse
import csv
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches

# Input rules and filtering
GROUP_SHEETS = None
ID_COL, EXPOSURE_COL, GENE_COL = "id", "Exposure", "Gene"
GROUP_ORDER = ["0h", "16h", "40h", "5d", "LA"]
UNSPECIFIED_GROUP = "Unspecified"
EXPOSURE_OVERRIDES = {}
EXCLUDED_RELATIONSHIPS = {("TN_NFKB", "USP14")}
GENE_ORDER_MODE = CELL_ORDER_MODE = "total_desc"
SHOW_ONLY_NONZERO_CELLS = True
SHOW_TOP_N_GENES = None

# Plot settings
DPI = 300
SAVE_PNG = SAVE_PDF = True
NODE_INCH = 0.2
DOT_RADIUS_MIN, DOT_RADIUS_MAX = 0.15, 0.48
BAR_WIDTH = 0.88
GROUP_COLORS = {
    "0h": "#6AA84F", "16h": "#F39C34", "40h": "#5B87B9",
    "5d": "#E76666", "LA": "#76A5AF", "Unspecified": "#999999",
}
FALLBACK_PALETTE = [
    "#59a14f", "#f28e2b", "#4e79a7", "#e15759", "#76b7b2",
    "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ab",
    "#c44764", "#2f9e6e", "#457b9d", "#b6992d", "#8cd17d",
]


def select_groups(observed, requested):
    groups = list(observed if requested is None else requested)
    if not groups or len(groups) != len(set(groups)):
        raise ValueError("Groups must be a nonempty list without duplicates.")
    missing = [g for g in groups if g not in observed]
    if missing:
        raise ValueError("Groups not found: " + ", ".join(missing))
    return groups


def load_csv(path, requested=None):
    for encoding in ("utf-8-sig", "gb18030"):
        try:
            with open(path, encoding=encoding, newline="") as handle:
                header = next((row for row in csv.reader(handle) if row), [])
            header = [name.strip() for name in header]
            if len(header) != len(set(header)):
                raise ValueError("Duplicate CSV column names after trimming spaces.")
            df = pd.read_csv(path, dtype=str, keep_default_na=False, encoding=encoding,
                             usecols=lambda name: name.strip() in (EXPOSURE_COL, GENE_COL))
            break
        except UnicodeDecodeError:
            if encoding == "gb18030":
                raise
    df.columns = df.columns.str.strip()
    required = [EXPOSURE_COL, GENE_COL]
    if not set(required).issubset(df.columns):
        raise ValueError("CSV requires columns: " + ", ".join(required))
    df = df[required].apply(lambda col: col.str.strip())
    df = df.loc[~df.eq("").all(axis=1)]
    incomplete = df.eq("").any(axis=1)
    if incomplete.any():
        raise ValueError("Missing Exposure/Gene in records: " + str((df.index[incomplete] + 2).tolist()[:10]))
    keep = [(exposure, gene) not in EXCLUDED_RELATIONSHIPS
            for exposure, gene in df.itertuples(index=False, name=None)]
    excluded = len(df) - sum(keep)
    df = df.loc[keep].copy()
    if df.empty:
        raise ValueError("No relationships remain after filtering.")
    mapping, unspecified = {}, []
    for exposure in df[EXPOSURE_COL].unique():
        if exposure in EXPOSURE_OVERRIDES:
            cell, group = map(str, EXPOSURE_OVERRIDES[exposure])
            cell, group = cell.strip(), group.strip()
        else:
            cell, separator, suffix = exposure.rpartition("_")
            if separator and suffix in GROUP_ORDER:
                group = suffix
            else:
                cell, group = exposure, UNSPECIFIED_GROUP
                unspecified.append(exposure)
        if not cell or not group:
            raise ValueError("Invalid cell/group for exposure: " + exposure)
        mapping[exposure] = (cell, group)
    df["Cell"] = df[EXPOSURE_COL].map(lambda value: mapping[value][0])
    df["Group"] = df[EXPOSURE_COL].map(lambda value: mapping[value][1])
    df = df.drop_duplicates([GENE_COL, "Cell", "Group"])
    observed = df["Group"].drop_duplicates().tolist()
    ordered = [g for g in GROUP_ORDER if g in observed]
    ordered += [g for g in observed if g not in ordered]
    groups = select_groups(ordered, requested)
    df = df.loc[df["Group"].isin(groups)]
    genes, cells = df[GENE_COL].unique().tolist(), df["Cell"].unique().tolist()
    mats = {}
    for group in groups:
        records = df.loc[df["Group"] == group]
        mats[group] = pd.crosstab(records[GENE_COL], records["Cell"]).reindex(
            index=genes, columns=cells, fill_value=0).clip(upper=1).astype(int)
    if unspecified:
        warnings.warn("No recognized time suffix; retained as " + UNSPECIFIED_GROUP
                      + ": " + ", ".join(unspecified), UserWarning)
    print("Excluded rows:", excluded)
    return groups, genes, cells, mats


def load_excel(path, requested=None):
    mats, genes, cells = {}, [], []
    with pd.ExcelFile(path) as workbook:
        available = workbook.sheet_names
        if requested is None:
            requested = [g for g in available if g.strip().lower() not in
                         {"merged_all", "merged", "merge", "all", "summary"}]
        groups = select_groups(available, requested)
        for group in groups:
            df = workbook.parse(group)
            df.columns = [str(c).strip() for c in df.columns]
            if df.columns.duplicated().any() or ID_COL not in df.columns:
                raise ValueError("Invalid or missing ID column in sheet: " + group)
            df[ID_COL] = df[ID_COL].fillna("").astype(str).str.strip()
            df = df.loc[df[ID_COL].ne("") & df[ID_COL].str.lower().ne("nan")]
            genes = list(dict.fromkeys(genes + df[ID_COL].tolist()))
            cell_cols = [c for c in df.columns if c != ID_COL]
            cells = list(dict.fromkeys(cells + cell_cols))
            values = df[cell_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
            values = values.astype(int).clip(0, 1)
            values.index = df[ID_COL]
            mats[group] = values.groupby(level=0, sort=False).max()
    mats = {g: mats[g].reindex(index=genes, columns=cells, fill_value=0) for g in groups}
    return groups, genes, cells, mats


def prepare_data(path, top_n=None):
    extension = Path(path).suffix.lower()
    if extension == ".csv":
        groups, genes, cells, mats = load_csv(path, GROUP_SHEETS)
    elif extension in {".xlsx", ".xls", ".xlsm"}:
        groups, genes, cells, mats = load_excel(path, GROUP_SHEETS)
    else:
        raise ValueError("Input must be CSV or Excel (.xlsx, .xls, .xlsm).")
    data = np.stack([mats[g].to_numpy(dtype=int) for g in groups])
    if top_n is not None:
        if isinstance(top_n, bool) or not isinstance(top_n, int) or top_n < 1:
            raise ValueError("Top N must be a positive integer.")
        top = np.argsort(-data.sum(axis=(0, 2)), kind="stable")[:top_n]
        keep = sorted(top)
        genes, data = [genes[i] for i in keep], data[:, keep, :]
    if SHOW_ONLY_NONZERO_CELLS:
        keep = np.flatnonzero(data.sum(axis=(0, 1)) > 0)
        cells, data = [cells[i] for i in keep], data[:, :, keep]
    if not genes or not cells:
        raise ValueError("No genes/cells to plot after filtering.")

    def order(items, totals, mode):
        if mode == "appearance":
            return np.arange(len(items))
        if mode == "alphabetical":
            return np.argsort(items, kind="stable")
        if mode == "total_desc":
            return np.argsort(-totals, kind="stable")
        raise ValueError("Unknown order mode: " + mode)

    gi = order(genes, data.sum(axis=(0, 2)), GENE_ORDER_MODE)
    ci = order(cells, data.sum(axis=(0, 1)), CELL_ORDER_MODE)
    return groups, [genes[i] for i in gi], [cells[i] for i in ci], data[:, gi, :][:, :, ci]


def set_count_ticks(ax, axis, end):
    end = max(1, int(end))
    getattr(ax, "set_" + axis + "lim")(0, end)
    raw = end / 4
    base = 10 ** np.floor(np.log10(raw))
    step = max(1, next(m * base for m in (1, 2, 5, 10) if m * base >= raw))
    ticks = sorted(set(np.arange(0, end + 1e-9, step).tolist() + [end]))
    labels = [str(int(t)) for t in ticks]
    extent = ax.get_window_extent()
    pixels = extent.width if axis == "x" else extent.height
    if len(ticks) > 2 and pixels * (ticks[-1] - ticks[-2]) / end < 45:
        labels[-2] = ""
    getattr(ax, "set_" + axis + "ticks")(ticks)
    getattr(ax, "set_" + axis + "ticklabels")(labels)


def plot_relationships(groups, genes, cells, data):
    count = len(groups)
    colors = [GROUP_COLORS.get(g, FALLBACK_PALETTE[i % len(FALLBACK_PALETTE)])
              for i, g in enumerate(groups)]
    radii = np.linspace(DOT_RADIUS_MIN, DOT_RADIUS_MAX, count) if count > 1 else [DOT_RADIUS_MAX]
    diameters = 2 * NODE_INCH * np.asarray(radii)
    ncols = 1 if count <= 4 else (2 if count <= 10 else 3)
    rows = int(np.ceil(count / ncols))
    col_width = max(1.1, max(len(g) for g in groups) * 0.085 + 0.35)
    legend_w = max(2.4, ncols * col_width + 0.2)
    size_rows = np.maximum(diameters, 0.18) + 0.08
    legend_min_h = 0.8 + rows * 0.25 + size_rows.sum()
    matrix_w, matrix_h = len(genes) * NODE_INCH, len(cells) * NODE_INCH
    top_h, right_w, gap = 4.6 * NODE_INCH, 8.5 * NODE_INCH, 0.06
    legend_h = max(matrix_h + gap + top_h, legend_min_h)
    left, bottom = 1.35, 1.35
    fig_w = left + matrix_w + 2 * gap + right_w + legend_w + 0.25
    fig_h = bottom + legend_h + 0.35
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=DPI)

    def axes(x, y, w, h):
        return fig.add_axes([x / fig_w, y / fig_h, w / fig_w, h / fig_h])

    matrix = axes(left, bottom, matrix_w, matrix_h)
    matrix.set(xlim=(-0.5, len(genes) - 0.5), ylim=(len(cells) - 0.5, -0.5), aspect="equal")
    matrix.set_xticks(range(len(genes)))
    matrix.set_xticklabels(genes, rotation=45, ha="right", fontsize=7)
    matrix.set_yticks(range(len(cells)))
    matrix.set_yticklabels(cells, fontsize=9)
    matrix.tick_params(axis="both", length=0)
    for spine in matrix.spines.values():
        spine.set_linewidth(1.8)
    for xi, yi in np.argwhere(data.any(axis=0)):
        present = np.flatnonzero(data[:, xi, yi])
        for k, group_index in enumerate(present):
            matrix.add_patch(patches.Wedge(
                (xi, yi), radii[len(present) - 1],
                90 + 360 * k / len(present), 90 + 360 * (k + 1) / len(present),
                facecolor=colors[group_index], edgecolor="white", linewidth=0.5))

    top = axes(left, bottom + matrix_h + gap, matrix_w, top_h)
    right = axes(left + matrix_w + gap, bottom, right_w, matrix_h)
    for ax, counts, horizontal in ((top, data.sum(axis=2), False), (right, data.sum(axis=1), True)):
        totals = np.zeros(counts.shape[1])
        positions = np.arange(counts.shape[1])
        for values, color in zip(counts, colors):
            if horizontal:
                ax.barh(positions, values, left=totals, height=BAR_WIDTH, color=color)
            else:
                ax.bar(positions, values, bottom=totals, width=BAR_WIDTH, color=color)
            totals += values
        ax.margins(0)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        if horizontal:
            ax.set_ylim(matrix.get_ylim())
            ax.set_yticks([])
            ax.set_xlabel("Count", fontsize=9)
        else:
            ax.set_xlim(matrix.get_xlim())
            ax.set_xticks([])
            ax.set_ylabel("Count", fontsize=9)
        set_count_ticks(ax, "x" if horizontal else "y", totals.max())

    # Legend coordinates are inches, so circle sizes match the matrix exactly.
    legend = axes(left + matrix_w + right_w + 2 * gap, bottom, legend_w, legend_h)
    legend.set(xlim=(0, legend_w), ylim=(0, legend_h))
    legend.axis("off")
    legend.text(0.08, legend_h - 0.08, "Group", fontsize=10, va="top")
    for i, (group, color) in enumerate(zip(groups, colors)):
        column, row = divmod(i, rows)
        x, y = 0.08 + column * col_width, legend_h - 0.4 - row * 0.25
        legend.add_patch(patches.Rectangle((x, y - 0.065), 0.16, 0.13, facecolor=color))
        legend.text(x + 0.24, y, group, fontsize=9, va="center")
    y = legend_h - 0.45 - rows * 0.25
    legend.text(0.08, y, "Number of groups", fontsize=10, va="top")
    y -= 0.25
    for level, diameter, height in zip(range(1, count + 1), diameters, size_rows):
        cy = y - height / 2
        legend.add_patch(patches.Circle((0.25, cy), diameter / 2, facecolor="#333333"))
        legend.text(0.52, cy, str(level), fontsize=9, va="center")
        y -= height
    return fig


def main(input_path, output_prefix=None, top_n=SHOW_TOP_N_GENES):
    groups, genes, cells, data = prepare_data(input_path, top_n)
    prefix = Path(output_prefix) if output_prefix else Path(__file__).with_name("upset_like_csv")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    print("Loaded: {} genes, {} cells, {} groups, {} relationships.".format(
        len(genes), len(cells), len(groups), int(data.sum())))
    fig = plot_relationships(groups, genes, cells, data)
    try:
        for extension, enabled in (("png", SAVE_PNG), ("pdf", SAVE_PDF)):
            if enabled:
                path = str(prefix) + "." + extension
                fig.savefig(path, bbox_inches="tight", pad_inches=0.02)
                print("Saved:", path)
    finally:
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="CSV/Excel input path")
    parser.add_argument("--output-prefix", help="Output path without extension")
    parser.add_argument("--top-n", type=int, default=SHOW_TOP_N_GENES, help="Show the top N genes")
    args = parser.parse_args()
    main(args.input, args.output_prefix, args.top_n)
