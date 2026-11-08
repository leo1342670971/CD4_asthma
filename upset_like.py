
import os
import csv
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.ticker import FuncFormatter


# ===================== USER SETTINGS（用户参数设置） =====================
# 1) 输入/输出路径
INPUT_PATH = r"C:\Users\13426\Desktop\工作簿2.csv"  # CSV: Exposure,Gene；也兼容原来的多sheet Excel
OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))             # 默认输出到本脚本所在文件夹
OUT_PNG    = os.path.join(OUTPUT_DIR, "upset_like_csv.png")
OUT_PDF    = os.path.join(OUTPUT_DIR, "upset_like_csv.pdf")
SAVE_PNG   = True                                                  # 是否保存PNG
SAVE_PDF   = True                                                  # 是否保存PDF
DPI = 300                                                          # 输出分辨率

# 2) Sheet / 数据结构设置
GROUP_SHEETS = None                                                # None=自动读取所有sheet；也可指定列表来控制顺序
ID_COL = "id"                                                      # 基因列名（每个sheet中该列存放基因名/ID）
EXPOSURE_COL = "Exposure"                                          # CSV中的细胞/时间组合列
GENE_COL = "Gene"                                                  # CSV中的基因列
GROUP_ORDER = ["0h", "16h", "40h", "5d", "LA"]                      # 识别的时间后缀及显示顺序
UNSPECIFIED_GROUP = "Unspecified"                                  # 缺少时间后缀的记录单独显示，不推测时间
EXPOSURE_OVERRIDES = {}                                            # 例如 {"TN_NFKB": ("TN_NFKB", "5d")}；仅确认后填写
EXCLUDED_RELATIONSHIPS = {("TN_NFKB", "USP14")}                      # 用户确认排除的(Exposure, Gene)关系；设为set()可恢复

# 3) 排序与筛选
GENE_ORDER_MODE = "total_desc"                                     # 基因排序：appearance=按首次出现顺序；alphabetical=字母序；total_desc=按总关系数降序
CELL_ORDER_MODE = "total_desc"                                     # 细胞排序：同上
SHOW_ONLY_NONZERO_CELLS = True                                     # True=只显示至少在某个group里出现过关系（计数>0）的细胞
SHOW_TOP_N_GENES = None                                          # 只显示前多少名的基因（按总关系数的多少计算），None=全部显示

# 4) 颜色（group配色）
GROUP_COLORS = {
    "0h":  "#6AA84F",
    "16h": "#F39C34",
    "40h": "#5B87B9",
    "5d":  "#E76666",
    "LA":  "#76A5AF",
    "Unspecified": "#999999"
}
FALLBACK_PALETTE = [                                               # 当sheet名不在GROUP_COLORS里时自动按顺序补色
    "#59a14f", "#f28e2b", "#4e79a7", "#e15759", "#76b7b2",
    "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ab",
    "#c44764", "#2f9e6e", "#457b9d", "#b6992d", "#8cd17d"
]

# 5) 画布“自适应布局”参数（通过节点间距自动推导画布大小、面板大小）
NODE_INCH = 0.2                                                   # 一个“节点格子”的物理尺寸（英寸），越大=节点间距越大=图越大
PAD_Y_INCH = 0.06                                                  # 上方直方图与矩阵之间的间隔
PAD_X_INCH = 0.06                                                  # 右侧直方图与矩阵之间的间隔
TOP_PANEL_NODES   = 4.6                                            # 上方直方图面板高度≈多少个节点格子
RIGHT_PANEL_NODES = 8.5                                            # 右侧直方图面板宽度≈多少个节点格子
LEGEND_PANEL_BASE_INCH = 3.6                                       # 右侧图例面板基础宽度（会随group数与文字长度自动扩展）
LEFT_MARGIN_INCH   = 1.35                                          # 左边留白（用于y轴标签）
BOTTOM_MARGIN_INCH = 1.35                                          # 底部留白（用于x轴标签）
TOP_MARGIN_INCH    = 0.35                                          # 顶部留白
RIGHT_MARGIN_INCH  = 0.25                                          # 最右留白

# 6) 直方图条宽/条高与条间隙（与节点格子对齐并留空隙）
BAR_GAP = 0.12                                                     # 条形间的空隙比例（相对一个格子宽/高）；越大=空隙越大
TOP_BAR_WIDTH     = 1.0 - BAR_GAP                                  # 上方每个基因柱的宽度（单位：节点格子）
RIGHT_BAR_HEIGHT  = 1.0 - BAR_GAP                                  # 右侧每个细胞条的高度（单位：节点格子）

# 7) 中间矩阵拼色圆点大小（单位：节点格子的“数据坐标半径”）
DOT_RADIUS_MIN = 0.15                                              # 只有1个group命中时的圆点半径
DOT_RADIUS_MAX = 0.48                                              # 多个group命中时的圆点最大半径（随命中group数线性插值）
WEDGE_EDGE_COLOR = "white"                                         # 拼色扇形之间的分隔线颜色
WEDGE_EDGE_LW    = 0.50                                            # 拼色扇形之间的分隔线宽度

# 8) 字体与标签
plt.rcParams["axes.unicode_minus"] = False                         # 防止负号显示异常
X_LABEL_ROT  = 45                                                  # 基因名倾斜角度
X_LABEL_SIZE = 7                                                   # 基因名字体大小
Y_LABEL_SIZE = 9                                                   # 细胞名字体大小
TOP_VALUE_AXIS_LABEL   = "Count"                                   # 上方直方图y轴标题
RIGHT_VALUE_AXIS_LABEL = "Count"                                   # 右侧直方图x轴标题

# 9) 线条粗细与颜色
MATRIX_BORDER_LW = 1.8                                             # 中间矩阵边框粗细
AXIS_LW = 1.0                                                      # 直方图坐标轴线粗细
BORDER_COLOR = "black"                                             # 边框/轴线颜色

# 10) 图例字体
LEGEND_FONT_SIZE  = 9                                              # 图例文字字号
LEGEND_TITLE_SIZE = 10                                             # 图例标题字号
SIZE_LEGEND_COLOR = "#333333"                                      # “Number of groups”圆的颜色

# 11) 刻度重叠处理（你最后提的需求）
TICK_END_MIN_PX = 45                                               # 当最高刻度(end)与相邻刻度在屏幕上太近时：保留刻度线，但隐藏相邻刻度数字
# ==========================================

def _safe_makedirs_for_file(fp: str) -> None:
    d = os.path.dirname(os.path.abspath(fp))
    if d and not os.path.exists(d):
        os.makedirs(d, exist_ok=True)


def _default_sheet_filter(name: str) -> bool:
    s = str(name).strip().lower()
    return s not in ("merged_all", "merged", "merge", "all", "summary")


def load_group_matrices(excel_path: str, group_sheets=None, id_col: str = "id"):
    xl = pd.ExcelFile(excel_path)
    all_sheets = xl.sheet_names

    if group_sheets is None:
        group_sheets = [s for s in all_sheets if _default_sheet_filter(s)]
    if len(group_sheets) < 1:
        raise ValueError("No group sheets found. Please set GROUP_SHEETS explicitly.")

    mats = {}
    gene_order, cell_order = [], []
    gene_seen, cell_seen = set(), set()

    for g in group_sheets:
        df = xl.parse(g)
        # 统一清理列名，避免 Excel 表头前后有空格导致 reindex 后整列变成 0
        df.columns = [str(c).strip() for c in df.columns]

        if id_col not in df.columns:
            raise ValueError(f"Sheet '{g}' missing required column '{id_col}'")

        df = df.copy()
        # 避免空 ID 被 astype(str) 转成字符串 "nan" 后作为一个基因进入作图
        df[id_col] = df[id_col].where(df[id_col].notna(), "").astype(str).str.strip()
        df = df[(df[id_col] != "") & (df[id_col].str.lower() != "nan")]

        for gene in df[id_col].tolist():
            if gene and gene not in gene_seen:
                gene_order.append(gene)
                gene_seen.add(gene)

        for c in df.columns:
            if c == id_col:
                continue
            c = str(c).strip()
            if c and c not in cell_seen:
                cell_order.append(c)
                cell_seen.add(c)

        cell_cols = [c for c in df.columns if c != id_col]
        for c in cell_cols:
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int).clip(0, 1)

        df = df.set_index(id_col)
        if df.index.duplicated().any():
            df = df.groupby(level=0).max()

        mats[g] = df

    genes = gene_order
    cells = cell_order

    for g in group_sheets:
        mats[g] = mats[g].reindex(index=genes, columns=cells, fill_value=0).astype(int)

    return group_sheets, genes, cells, mats


def load_csv_matrices(csv_path: str, group_sheets=None):
   
    frames = []
    blank_rows = 0
    for encoding in ("utf-8-sig", "gb18030"):
        frames = []
        blank_rows = 0
        try:
            # pandas会自动重命名重复表头，因此在读取前检查原始表头。
            with open(csv_path, "r", encoding=encoding, newline="") as handle:
                header = next((row for row in csv.reader(handle) if row), [])
                header = [name.strip() for name in header]
                if len(header) != len(set(header)):
                    raise ValueError("CSV has duplicate column names after trimming spaces.")
            with pd.read_csv(csv_path, dtype=str, keep_default_na=False,
                             encoding=encoding, chunksize=100000) as reader:
                for chunk in reader:
                    chunk.columns = [str(c).strip() for c in chunk.columns]
                    required = [EXPOSURE_COL, GENE_COL]
                    missing = [c for c in required if c not in chunk.columns]
                    if missing:
                        raise ValueError("CSV missing required columns: " + ", ".join(missing))
                    if chunk.columns.duplicated().any():
                        raise ValueError("CSV has duplicate column names after trimming spaces.")
                    chunk = chunk.loc[:, required].copy()
                    for col in required:
                        chunk[col] = chunk[col].str.strip()
                    blank = (chunk == "").all(axis=1)
                    blank_rows += int(blank.sum())
                    chunk = chunk.loc[~blank]
                    incomplete = (chunk == "").any(axis=1)
                    if incomplete.any():
                        rows = (chunk.index[incomplete] + 2).tolist()[:10]
                        raise ValueError(f"CSV rows missing Exposure or Gene (header=1): {rows}")
                    if not chunk.empty:
                        frames.append(chunk)
            break
        except UnicodeDecodeError:
            if encoding == "gb18030":
                raise
    if not frames:
        raise ValueError("CSV contains no valid Exposure/Gene records.")

    df = pd.concat(frames, ignore_index=True)
    row_count = len(df)
    excluded = pd.Series(
        [(exposure, gene) in EXCLUDED_RELATIONSHIPS
         for exposure, gene in df[[EXPOSURE_COL, GENE_COL]].itertuples(index=False, name=None)],
        index=df.index,
    )
    excluded_rows = int(excluded.sum())
    df = df.loc[~excluded].copy()
    if df.empty:
        raise ValueError("No relationships remain after applying EXCLUDED_RELATIONSHIPS.")
    mapping, unlabelled = {}, []
    for exposure in df[EXPOSURE_COL].unique():
        if exposure in EXPOSURE_OVERRIDES:
            cell, group = EXPOSURE_OVERRIDES[exposure]
            cell, group = str(cell).strip(), str(group).strip()
        else:
            cell, separator, suffix = exposure.rpartition("_")
            if separator and suffix in GROUP_ORDER:
                group = suffix
            else:
                cell, group = exposure, UNSPECIFIED_GROUP
                unlabelled.append(exposure)
        if not cell or not group:
            raise ValueError(f"Invalid cell/group for Exposure: {exposure!r}")
        mapping[exposure] = (cell, group)
    df["Cell"] = df[EXPOSURE_COL].map(lambda value: mapping[value][0])
    df["Group"] = df[EXPOSURE_COL].map(lambda value: mapping[value][1])
    df = df.drop_duplicates(subset=[GENE_COL, "Cell", "Group"])
    duplicates = row_count - excluded_rows - len(df)

    observed_groups = df["Group"].drop_duplicates().tolist()
    if group_sheets is None:
        groups = [g for g in GROUP_ORDER if g in observed_groups]
        groups += [g for g in observed_groups if g not in groups]
    else:
        groups = list(group_sheets)
        if not groups or len(groups) != len(set(groups)):
            raise ValueError("GROUP_SHEETS must be a nonempty list without duplicates.")
        missing = [g for g in groups if g not in observed_groups]
        if missing:
            raise ValueError(f"Groups not found in CSV: {missing}")
        df = df.loc[df["Group"].isin(groups)]
    genes = df[GENE_COL].drop_duplicates().tolist()
    cells = df["Cell"].drop_duplicates().tolist()
    mats = {}
    for group in groups:
        records = df.loc[df["Group"] == group]
        mat = pd.DataFrame(0, index=genes, columns=cells, dtype=int)
        for gene, cell in records[[GENE_COL, "Cell"]].itertuples(index=False, name=None):
            mat.at[gene, cell] = 1
        mats[group] = mat

    if unlabelled:
        warnings.warn("No recognized time suffix; retained in '" + UNSPECIFIED_GROUP
                      + "': " + ", ".join(unlabelled), UserWarning)
    print(f"CSV: {row_count} nonblank rows; {blank_rows} blank rows skipped; "
          f"{excluded_rows} rows excluded; {duplicates} duplicate relationships merged.")
    print(f"Loaded: {len(genes)} genes, {len(cells)} cells, {len(groups)} groups, "
          f"{len(df)} unique relationships.")
    print("Group relationship counts:", {g: int(mats[g].to_numpy().sum()) for g in groups})
    return groups, genes, cells, mats


def load_input_matrices(input_path: str, group_sheets=None, id_col: str = "id"):
    extension = os.path.splitext(os.fspath(input_path))[1].lower()
    if extension == ".csv":
        return load_csv_matrices(input_path, group_sheets)
    if extension in (".xlsx", ".xls", ".xlsm"):
        return load_group_matrices(input_path, group_sheets, id_col)
    raise ValueError("Input must be a .csv, .xlsx, .xls or .xlsm file.")


def reorder_items(order_mode: str, items: list, totals: pd.Series):
    if order_mode == "appearance":
        return items
    if order_mode == "alphabetical":
        return sorted(items)
    if order_mode == "total_desc":
        return totals.loc[items].sort_values(ascending=False, kind="stable").index.tolist()
    raise ValueError(f"Unknown order mode: {order_mode}")


def _int_formatter():
    return FuncFormatter(lambda x, pos: f"{int(x)}" if abs(x - int(x)) < 1e-6 else f"{x:g}")


def _nice_ceil_int(v: float) -> int:
    return int(np.ceil(float(v)))


def _nice_step(raw_step: float) -> float:
    if raw_step <= 0:
        return 1.0
    k = np.floor(np.log10(raw_step))
    base = 10 ** k
    for m in (1, 2, 5, 10):
        step = m * base
        if step >= raw_step:
            return step
    return 10 * base


def _set_nice_ticks_0_to_end(ax, axis: str, end_value: float, nbins: int = 5, min_px: int = 14):
    end_value = float(end_value)
    if end_value <= 0:
        if axis == "x":
            ax.set_xticks([0])
            ax.set_xticklabels(["0"])
        else:
            ax.set_yticks([0])
            ax.set_yticklabels(["0"])
        return

    intervals = max(1, nbins - 1)
    raw_step = end_value / intervals
    step = max(1.0, _nice_step(raw_step))  # 关系数为整数，避免出现0.2/0.5刻度

    ticks = list(np.arange(0, end_value + 1e-9, step))
    if abs(ticks[-1] - end_value) > 1e-6:
        ticks.append(end_value)

    ticks = [t for t in ticks if t <= end_value + 1e-6]
    if not any(abs(t - end_value) < 1e-6 for t in ticks):
        ticks.append(end_value)
    ticks = sorted(set(ticks))

    hide_labels = set()

    ax.figure.canvas.draw()
    trans = ax.transData
    if len(ticks) >= 2:
        prev = ticks[-2]
        if axis == "x":
            px_end = trans.transform((end_value, 0))[0]
            px_prev = trans.transform((prev, 0))[0]
        else:
            px_end = trans.transform((0, end_value))[1]
            px_prev = trans.transform((0, prev))[1]
        if abs(px_end - px_prev) < min_px:
            hide_labels.add(prev)

    fmt = _int_formatter()
    labels = []
    for t in ticks:
        if any(abs(t - h) < 1e-6 for h in hide_labels):
            labels.append("")
        else:
            labels.append(fmt(t, None))

    if axis == "x":
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
    else:
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)


def _pick_group_legend_ncol(G: int) -> int:
    if G <= 4:
        return 1
    if G <= 10:
        return 2
    return 3


def _radius_data_to_diameter_px(ax, r_data: float) -> float:
    fig = ax.figure
    fig.canvas.draw()
    bbox = ax.get_window_extent(fig.canvas.get_renderer())
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    px_per_data_x = bbox.width / max(1e-9, (xlim[1] - xlim[0]))
    px_per_data_y = bbox.height / max(1e-9, abs(ylim[1] - ylim[0]))
    px_per_data = min(px_per_data_x, px_per_data_y)
    return float(2.0 * r_data * px_per_data)


def _draw_size_legend_manual(fig, ax_leg, title, levels, diam_px_list, y_top):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axbb = ax_leg.get_window_extent(renderer)
    ax_h_px = axbb.height

    diam_frac = [d / ax_h_px for d in diam_px_list]

    title_h = 0.055
    edge_gap = 0.020
    need = title_h + sum(diam_frac) + edge_gap * (len(levels) - 1)

    avail = max(0.10, y_top)
    scale = 1.0
    if need > avail:
        free = max(0.05, avail - title_h)
        base = sum(diam_frac) + edge_gap * (len(levels) - 1)
        if base > 1e-9:
            scale = min(1.0, free / base)

    diam_frac = [d * scale for d in diam_frac]
    edge_gap *= scale

    x_circle = 0.12
    x_text   = 0.22

    ax_leg.text(0.08, y_top, title, ha="left", va="top",
                fontsize=LEGEND_TITLE_SIZE, transform=ax_leg.transAxes)

    y = y_top - title_h
    for lvl, d in zip(levels, diam_frac):
        r = d / 2.0
        cy = y - r
        ax_leg.add_patch(
            patches.Ellipse((x_circle, cy), width=d * ax_h_px / axbb.width, height=d,
                           transform=ax_leg.transAxes,
                           facecolor=SIZE_LEGEND_COLOR, edgecolor="none")
        )
        ax_leg.text(x_text, cy, str(lvl), ha="left", va="center",
                    fontsize=LEGEND_FONT_SIZE, transform=ax_leg.transAxes)
        y = y - d - edge_gap


def main(input_path=None, output_prefix=None):
    input_path = INPUT_PATH if input_path is None else input_path
    out_png = OUT_PNG if output_prefix is None else os.fspath(output_prefix) + ".png"
    out_pdf = OUT_PDF if output_prefix is None else os.fspath(output_prefix) + ".pdf"
    groups, genes, cells, mats = load_input_matrices(input_path, GROUP_SHEETS, ID_COL)
    G = len(groups)

    group_colors = dict(GROUP_COLORS)
    for i, g in enumerate(groups):
        if g not in group_colors:
            group_colors[g] = FALLBACK_PALETTE[i % len(FALLBACK_PALETTE)]

    gene_counts = pd.DataFrame({g: mats[g].sum(axis=1) for g in groups}, index=genes)
    cell_counts = pd.DataFrame({g: mats[g].sum(axis=0) for g in groups}, index=cells)

    if SHOW_TOP_N_GENES is not None and (not isinstance(SHOW_TOP_N_GENES, int) or SHOW_TOP_N_GENES < 1):
        raise ValueError("SHOW_TOP_N_GENES must be a positive integer or None.")
    if SHOW_TOP_N_GENES is not None and SHOW_TOP_N_GENES < len(genes):
        total_gene = gene_counts.sum(axis=1)
        top_set = set(total_gene.sort_values(ascending=False, kind="stable").head(SHOW_TOP_N_GENES).index)
        genes = [gn for gn in genes if gn in top_set]
        gene_counts = gene_counts.loc[genes]
        for g in groups:
            mats[g] = mats[g].loc[genes, :]
        # 如果只显示 Top N 基因，右侧细胞柱状图也应基于过滤后的基因重新计数
        cell_counts = pd.DataFrame({g: mats[g].sum(axis=0) for g in groups}, index=cells)

    if SHOW_ONLY_NONZERO_CELLS:
        keep_cells = cell_counts.sum(axis=1) > 0
        cells = [c for c in cells if bool(keep_cells.loc[c])]
        cell_counts = cell_counts.loc[cells]
        for g in groups:
            mats[g] = mats[g].loc[:, cells]

    genes = reorder_items(GENE_ORDER_MODE, genes, gene_counts.sum(axis=1))
    cells = reorder_items(CELL_ORDER_MODE, cells, cell_counts.sum(axis=1))

    gene_counts = gene_counts.loc[genes]
    cell_counts = cell_counts.loc[cells]
    for g in groups:
        mats[g] = mats[g].loc[genes, cells]

    n_genes = len(genes)
    n_cells = len(cells)

    if n_genes == 0:
        raise ValueError("No genes to plot after filtering.")
    if n_cells == 0:
        raise ValueError("No cells to plot after filtering. Please check whether the input matrix contains any 1 values.")

    def count_to_radius(k: int) -> float:
        if G <= 1:
            return DOT_RADIUS_MAX
        k = max(1, min(k, G))
        return DOT_RADIUS_MIN + (DOT_RADIUS_MAX - DOT_RADIUS_MIN) * (k - 1) / (G - 1)

    max_label_len = max(len(str(g)) for g in groups) if groups else 10
    legend_w = max(LEGEND_PANEL_BASE_INCH, 2.6 + 0.22 * G + 0.06 * max_label_len)

    matrix_w = n_genes * NODE_INCH
    matrix_h = n_cells * NODE_INCH
    top_h    = TOP_PANEL_NODES * NODE_INCH
    right_w  = RIGHT_PANEL_NODES * NODE_INCH

    fig_w = LEFT_MARGIN_INCH + matrix_w + PAD_X_INCH + right_w + PAD_X_INCH + legend_w + RIGHT_MARGIN_INCH
    fig_h = BOTTOM_MARGIN_INCH + matrix_h + PAD_Y_INCH + top_h + TOP_MARGIN_INCH
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=DPI)

    def frac(x_inch, y_inch, w_inch, h_inch):
        return [x_inch / fig_w, y_inch / fig_h, w_inch / fig_w, h_inch / fig_h]

    mat_left = LEFT_MARGIN_INCH
    mat_bot  = BOTTOM_MARGIN_INCH
    top_left = mat_left
    top_bot  = mat_bot + matrix_h + PAD_Y_INCH
    right_left = mat_left + matrix_w + PAD_X_INCH
    right_bot  = mat_bot
    leg_left = right_left + right_w + PAD_X_INCH
    leg_bot  = mat_bot
    leg_h    = matrix_h + PAD_Y_INCH + top_h

    x = np.arange(n_genes)
    y = np.arange(n_cells)

    ax_mat = fig.add_axes(frac(mat_left, mat_bot, matrix_w, matrix_h))
    ax_mat.set_xlim(-0.5, n_genes - 0.5)
    ax_mat.set_ylim(-0.5, n_cells - 0.5)
    ax_mat.invert_yaxis()
    ax_mat.set_aspect("equal", adjustable="box")

    for sp in ax_mat.spines.values():
        sp.set_linewidth(MATRIX_BORDER_LW)
        sp.set_color(BORDER_COLOR)

    ax_mat.set_xticks(x)
    ax_mat.set_xticklabels(genes, rotation=X_LABEL_ROT, ha="right", va="top", fontsize=X_LABEL_SIZE)
    ax_mat.set_yticks(y)
    ax_mat.set_yticklabels(cells, fontsize=Y_LABEL_SIZE)
    ax_mat.tick_params(axis="both", length=0)

    for yi in range(n_cells):
        for xi in range(n_genes):
            present = [g for g in groups if mats[g].iat[xi, yi] == 1]
            if not present:
                continue
            r = count_to_radius(len(present))
            step = 360.0 / len(present)
            start = 90.0
            for k, gname in enumerate(present):
                ax_mat.add_patch(
                    patches.Wedge(
                        center=(xi, yi),
                        r=r,
                        theta1=start + k * step,
                        theta2=start + (k + 1) * step,
                        facecolor=group_colors[gname],
                        edgecolor=WEDGE_EDGE_COLOR,
                        linewidth=WEDGE_EDGE_LW
                    )
                )

    ax_top = fig.add_axes(frac(top_left, top_bot, matrix_w, top_h))
    ax_top.set_xlim(ax_mat.get_xlim())
    ax_top.margins(x=0, y=0)

    bottom = np.zeros(n_genes)
    for g in groups:
        h = gene_counts[g].to_numpy()
        ax_top.bar(x, h, bottom=bottom, color=group_colors[g], width=TOP_BAR_WIDTH, edgecolor="none")
        bottom += h

    ax_top.tick_params(axis="x", labelbottom=False, bottom=False)
    ax_top.set_ylabel(TOP_VALUE_AXIS_LABEL, fontsize=9)
    ax_top.yaxis.set_major_formatter(_int_formatter())
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.spines["left"].set_linewidth(AXIS_LW)
    ax_top.spines["bottom"].set_linewidth(AXIS_LW)
    ax_top.minorticks_off()

    vmax_top = float(gene_counts.sum(axis=1).max())
    endy = float(_nice_ceil_int(vmax_top))
    ax_top.set_ylim(0, endy)
    _set_nice_ticks_0_to_end(ax_top, "y", endy, nbins=5, min_px=TICK_END_MIN_PX)

    ax_right = fig.add_axes(frac(right_left, right_bot, right_w, matrix_h))
    ax_right.set_ylim(ax_mat.get_ylim())

    ax_right.margins(x=0, y=0)

    left = np.zeros(n_cells)
    for g in groups:
        w = cell_counts[g].to_numpy()
        ax_right.barh(y, w, left=left, color=group_colors[g], height=RIGHT_BAR_HEIGHT, edgecolor="none")
        left += w

    ax_right.tick_params(axis="y", which="both", left=False, right=False, labelleft=False)
    ax_right.set_yticks([])
    ax_right.set_xlabel(RIGHT_VALUE_AXIS_LABEL, fontsize=9)
    ax_right.xaxis.set_major_formatter(_int_formatter())
    ax_right.spines["top"].set_visible(False)
    ax_right.spines["right"].set_visible(False)
    ax_right.spines["left"].set_linewidth(AXIS_LW)
    ax_right.spines["bottom"].set_linewidth(AXIS_LW)
    ax_right.minorticks_off()

    vmax_right = float(cell_counts.sum(axis=1).max())
    endx = float(_nice_ceil_int(vmax_right))
    ax_right.set_xlim(0, endx)
    _set_nice_ticks_0_to_end(ax_right, "x", endx, nbins=5, min_px=TICK_END_MIN_PX)

    ax_leg = fig.add_axes(frac(leg_left, leg_bot, legend_w, leg_h))
    ax_leg.axis("off")

    ncol = _pick_group_legend_ncol(G)
    group_handles = [patches.Patch(facecolor=group_colors[g], edgecolor="none", label=g) for g in groups]
    leg1 = ax_leg.legend(
        handles=group_handles,
        title="Group",
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        ncol=ncol,
        frameon=False,
        fontsize=LEGEND_FONT_SIZE,
        title_fontsize=LEGEND_TITLE_SIZE,
        borderaxespad=0.0,
        labelspacing=0.6,
        handlelength=1.6,
        columnspacing=1.2
    )
    ax_leg.add_artist(leg1)

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bb1 = leg1.get_window_extent(renderer)
    axbb = ax_leg.get_window_extent(renderer)
    bb1_bottom_y_frac = (bb1.y0 - axbb.y0) / axbb.height
    y_top = max(0.10, bb1_bottom_y_frac - 0.05)

    levels = list(range(1, G + 1))
    diam_px = [_radius_data_to_diameter_px(ax_mat, count_to_radius(k)) for k in levels]
    _draw_size_legend_manual(fig, ax_leg, "Number of groups", levels, diam_px, y_top)

    if SAVE_PNG:
        _safe_makedirs_for_file(out_png)
        fig.savefig(out_png, bbox_inches="tight", pad_inches=0.02)
        print("Saved PNG:", out_png)
    if SAVE_PDF:
        _safe_makedirs_for_file(out_pdf)
        fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.02)
        print("Saved PDF:", out_pdf)

    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot gene/cell relationships from Exposure,Gene CSV or multi-sheet Excel.")
    parser.add_argument("input", nargs="?", default=INPUT_PATH, help="Input CSV/Excel path")
    parser.add_argument("--output-prefix", help="Output path without extension; saves PNG and PDF")
    args = parser.parse_args()
    main(args.input, args.output_prefix)
