# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 16:30:04 2025

@author: Vasiliki
"""
#%%
dataset_to_animal = {
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_1.vzg2": "mp16_1", #ko1a
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_3.vzg2": "mp16_1", #ko1b
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_0.vzg2": "mp16_1", #ko1c
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_5.vzg2": "wt_16_1", #ko2a
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_4.vzg2": "wt_16_1", #ko2b
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_2.vzg2": "wt_16_1", #ko2c
    
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_5.vzg2": "wt_16_2", #an1a
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_3.vzg2": "wt_16_2", #an1b
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_4.vzg2": "wt_16_2", #an1c
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_1.vzg2": "mp16_2", #an2a
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_2.vzg2": "mp16_2", #an2b
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_0.vzg2": "mp16_2", #an2c
  
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_4.vzg2": "mp11_1", #ko1
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_3.vzg2": "mp11_1", #ko2
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_0.vzg2": "mp11_1", #ko3
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_5.vzg2": "wt11_1", #wt1
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_1.vzg2": "wt11_1", #wt2
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_2.vzg2": "wt11_1", #wt3
 
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_1.vzg2": "wt11_2", #ko1
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_5.vzg2": "wt11_2", #ko2
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_0.vzg2": "wt11_2", #ko3
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_4.vzg2": "wt11_6", #wt1
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_3.vzg2": "wt11_6", #wt2
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_2.vzg2": "wt11_6", #wt3
    
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_3.vzg2": "mp11_3", #ko1
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_4.vzg2": "mp11_3", #ko2
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_0.vzg2": "wt11_3", #ko3
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_1.vzg2": "wt11_3", #wt1
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_2.vzg2": "wt11_3", #wt2
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_5.vzg2": "mp11_3", #wt3
    
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_0.vzg2": "mp8_2", #ko8c
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_1.vzg2": "mp8_1", #ko81
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_5.vzg2": "mp8_1", #ko82
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_2.vzg2": "mp8_1", #ko83
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_3.vzg2": "mp8_2", #ko8a
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_4.vzg2": "mp8_2", #ko8b    
        
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_1.vzg2": "mp8_3", #ko8a
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_5.vzg2": "mp8_3", #ko8b
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_0.vzg2": "mp8_3", #ko8c
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_3.vzg2": "mp16_4", #ko16a
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_2.vzg2": "mp16_4", #ko16b
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_4.vzg2": "mp16_4", #ko16c
    
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_5.vzg2": "mp11_4", #ko10a
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_3.vzg2": "mp11_4", #ko10b
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_4.vzg2": "mp11_4", #ko10c
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_2.vzg2": "mp16_3", #ko16a
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_1.vzg2": "mp16_3", #ko16b
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_0.vzg2": "mp16_3"  #ko16c
}

#%%
from __future__ import annotations
import os, re, json
from pathlib import Path
from typing import Optional, Tuple, List, Union
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import tifffile
import matplotlib.colors as mcolors

# ===================== USER TWEAKS =====================
SUBCLUSTER_H5AD: Optional[str] = r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad"

# Gene heatmap settings
GENE_NAME: Optional[str] = None         # e.g., "Th"; if None you'll be prompted
EXPR_SOURCE: str = "X"                  # "X", "raw", or a layer name (e.g., "log1p","counts","scaled")
LOG1P: bool = False                     # Apply log1p before plotting (if your source is raw counts)

# Robust color scaling (percentiles); set either to None to use min/max
VMIN_PERC: Optional[float] = 1.0
VMAX_PERC: Optional[float] = 99.0
CMAP = mcolors.LinearSegmentedColormap.from_list(
    "gray_to_bordeaux",
    ["#9c9c9c", "#BDBDBD", "#800020", "#2B000F"]  # low → high expression
)
# Region/label columns (will auto-detect if not found)
REGION_COL_HINTS = ["Datasets", "section_id", "Region", "region", "slide_id"]
LABEL_COL_HINTS  = ["leiden_sub_Th_Slc6a3", "leiden_sub", "DA_subclusters_Th_Slc6a3", "leiden"]

# Optional faint category overlay on top of heatmap
PLOT_CATEGORIES_OVER = False
LEGEND = True
MAX_LEGEND = 25

# DAPI / image
START_DAPI_DIR = r"Z:\Vasska\Merscope_data\202507021142_MITOPARK1_VMSC13902\region_wt1\images"
DAPI_FILE      = r"Z:\Vasska\Merscope_data\202507021142_MITOPARK1_VMSC13902\region_wt1\images\mosaic_DAPI_z3.tif"
PYRAMID_LEVEL = -1      # 0 = full-res, -1 = lowest level

# Plot look
DOT_SIZE = 8
ALPHA = 0.90
DPI = 180
OUT_DIR = r"F:/vizgen_clustering_out/overlays_one_section"

# Coordinates / transform
UM_PER_PIXEL_FALLBACK: Optional[float] = 0.108  # used only if coords look like microns and no transform is provided
INVERT_Y = False
MAX_POINTS = 400_000

# Transform control
FORCE_TRANSFORM_DIR: Optional[str] = "forward"  # "forward", "inverse", or None to auto
SKIP_TRANSFORM: bool = False               # True if xy are already mosaic pixels

# Color scale control
GLOBAL_COLOR_SCALE: bool = True   # one scale for all selected sections in this run
FIXED_VMIN: Optional[float] = 0.0  # set to a number to lock, e.g., 0.0
FIXED_VMAX: Optional[float] = 5.0  # set to a number to lock, e.g., 5.0



# ======================================================


# -------------- tiny GUI pickers --------------
def _pick_file(title: str, exts: Tuple[str,...]) -> Optional[Path]:
    try:
        from tkinter import Tk, filedialog
        root = Tk(); root.withdraw()
        fp = filedialog.askopenfilename(
            title=title,
            initialdir=START_DAPI_DIR if "DAPI" in title and START_DAPI_DIR else None,
            filetypes=[(f"*{e}", f"*{e}") for e in exts] + [("All files","*.*")]
        )
        root.destroy()
        return Path(fp) if fp else None
    except Exception:
        return None

def _pick_folder(title: str) -> Optional[Path]:
    try:
        from tkinter import Tk, filedialog
        root = Tk(); root.withdraw()
        fp = filedialog.askdirectory(title=title, initialdir=START_DAPI_DIR or None)
        root.destroy()
        return Path(fp) if fp else None
    except Exception:
        return None


# -------------- helpers --------------
def _find_xy_cols(df: pd.DataFrame) -> Tuple[str,str]:
    cands = [("global_x","global_y"),("x_global","y_global"),("x","y"),
             ("center_x","center_y"),("x_centroid","y_centroid")]
    lower = {c.lower(): c for c in df.columns}
    for a,b in cands:
        if a in lower and b in lower: return lower[a], lower[b]
    for a,b in [("X_obsm0","Y_obsm1")]:
        if a in df.columns and b in df.columns: return a,b
    raise KeyError("Could not find XY columns (looked for global_x/global_y, x_global/y_global, center_x/center_y, x/y).")

def _choose_regions(df: pd.DataFrame, region_col: str) -> List[str]:
    vals = sorted(pd.Series(df[region_col].astype(str).unique()).tolist())
    print("\nAvailable regions:")
    for i, v in enumerate(vals):
        print(f"{i:>3}  {v}")
    sel = input("\nType indices for regions (e.g., 0,2,5) or 'all' (default 0): ").strip()
    if not sel:
        return [vals[0]]
    if sel.lower() in ("all", "a", "*"):
        return vals
    try:
        idxs = [int(s) for s in re.split(r"[,\s;]+", sel) if s != ""]
        return [vals[i] for i in idxs]
    except Exception:
        toks = [s for s in re.split(r"[,\s;]+", sel) if s]
        picked = []
        for t in toks:
            hit = [v for v in vals if t.lower() in v.lower()]
            picked.extend(hit or [])
        return sorted(set(picked)) or [vals[0]]

def _read_pyramid_image(path: Path, pyramid_level: int = -1):
    with tifffile.TiffFile(str(path)) as tf:
        s = tf.series[0]
        if getattr(s, "levels", None) and len(s.levels) > 1:
            lvls = s.levels
            idx = pyramid_level if pyramid_level >= 0 else len(lvls) + pyramid_level
            idx = max(0, min(idx, len(lvls)-1))
            img = np.squeeze(lvls[idx].asarray())
            if img.ndim > 2: img = img[...,0]
            H_full, W_full = lvls[0].shape[0], lvls[0].shape[1]
            return img, (W_full, H_full), idx
        # single level → lightweight fallback
        try:
            import zarr
            arr = zarr.open(s.aszarr(), mode='r')
            H_full, W_full = int(arr.shape[0]), int(arr.shape[1])
            stride = max(1, int(round(max(W_full, H_full)/10000)))
            img = np.asarray(arr[::stride, ::stride])
            return img, (W_full, H_full), f"decim/{stride}"
        except Exception:
            page = np.squeeze(s.pages[0].asarray())
            if page.ndim > 2: page = page[...,0]
            H_full, W_full = page.shape[0], page.shape[1]
            return page, (W_full, H_full), 0

def _scale_to_level(x_full: np.ndarray, y_full: np.ndarray, W_full: int, H_full: int, W_lvl: int, H_lvl: int):
    sx, sy = W_lvl/float(W_full), H_lvl/float(H_full)
    return x_full*sx, y_full*sy

def _choose_palette(n: int):
    import matplotlib
    if n <= 10:  cmap = matplotlib.cm.get_cmap("tab10", n)
    elif n <= 20: cmap = matplotlib.cm.get_cmap("tab20", n)
    else:         cmap = matplotlib.cm.get_cmap("hsv", n)
    return [cmap(i) for i in range(n)]

def _augment_to_3x3(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, float)
    if M.shape == (2,3):
        M = np.vstack([M, [0,0,1]])
    if M.shape != (3,3):
        raise ValueError(f"Transform must be 2x3 or 3x3, got {M.shape}")
    return M

def _load_transform_matrix(path: Path) -> np.ndarray:
    ext = path.suffix.lower()
    if ext in (".json",".vzg",".vzg2"):
        data = json.load(open(path, "r", encoding="utf-8"))
        keys = ["micron_to_mosaic_pixel_transform","micron_to_mosaic_pixel_matrix",
                "microns_to_mosaic_pixels","mosaic_to_micron","mosaic_pixel_to_micron_transform"]
        def _find(obj):
            if isinstance(obj, dict):
                for k,v in obj.items():
                    if k in keys: return v
                    hit = _find(v)
                    if hit is not None: return hit
            if isinstance(obj, list):
                for v in obj:
                    hit = _find(v)
                    if hit is not None: return hit
            return None
        M = _find(data)
        if M is None:
            raise KeyError(f"No known transform key in {path.name}")
        return _augment_to_3x3(np.asarray(M, float))

    # delimiter-aware numeric reader for .txt/.csv (whitespace/comma/semicolon)
    try:
        arr = np.genfromtxt(str(path), delimiter=None)
        if not np.isfinite(arr).all():
            arr = np.genfromtxt(str(path), delimiter=",")
        if not np.isfinite(arr).all():
            arr = np.genfromtxt(str(path), delimiter=";")
    except Exception:
        arr = np.loadtxt(str(path), dtype=float)

    flat = np.asarray(arr, float).ravel()
    if flat.size == 6:
        M = flat.reshape(2,3)
    elif flat.size == 9:
        M = flat.reshape(3,3)
    else:
        raise ValueError(f"Transform file must contain 6 or 9 numbers, got {flat.size}.")
    return _augment_to_3x3(M)

def _apply_affine(M: np.ndarray, x_um: np.ndarray, y_um: np.ndarray) -> Tuple[np.ndarray,np.ndarray]:
    P = np.vstack([x_um, y_um, np.ones_like(x_um, float)])
    Q = (M @ P).T
    return Q[:,0], Q[:,1]

def _score_mapping(x_pix: np.ndarray, y_pix: np.ndarray, W: int, H: int):
    inside = (x_pix >= -0.05*W) & (x_pix <= 1.05*W) & (y_pix >= -0.05*H) & (y_pix <= 1.05*H)
    frac = float(np.mean(inside)) if inside.size else 0.0
    xs, ys = x_pix[inside], y_pix[inside]
    if xs.size > 0:
        spread = float(np.sqrt(np.var(xs) + np.var(ys)))
        cx, cy = float(np.nanmean(xs)), float(np.nanmean(ys))
        center_dist = float(np.hypot(cx - W/2.0, cy - H/2.0))
    else:
        spread, center_dist = 0.0, 1e12
    score = frac * (1.0 + spread / max(W, H)) - 0.1 * (center_dist / max(W, H))
    return score, frac, spread, center_dist

def _auto_transform_direction(M: np.ndarray, x_um: np.ndarray, y_um: np.ndarray, W: int, H: int):
    from numpy.linalg import inv
    xf, yf = _apply_affine(M, x_um, y_um)
    sf, fracf, spf, cdf = _score_mapping(xf, yf, W, H)
    try:
        Mi = inv(M)
        xi, yi = _apply_affine(Mi, x_um, y_um)
        si, fraci, spi, cdi = _score_mapping(xi, yi, W, H)
    except Exception:
        si, fraci, spi, cdi = -1e9, -1.0, 0.0, 1e12
        xi, yi = np.array([]), np.array([])

    sn, fracn, spn, cdn = _score_mapping(x_um, y_um, W, H)  # none
    scores = [
        ("forward", sf, xf, yf, fracf),
        ("inverse", si, xi, yi, fraci),
        ("none",    sn, x_um, y_um, fracn),
    ]
    best = max(scores, key=lambda t: t[1])
    return best  # (label, score, x, y, frac_inside)

def _infer_unit(x: np.ndarray, y: np.ndarray, W: int, H: int) -> str:
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if x.size == 0: return "pixel"
    xq, yq = np.quantile(np.abs(x), 0.95), np.quantile(np.abs(y), 0.95)
    return "pixel" if (xq <= W*1.5 and yq <= H*1.5) else "micron"


# -------------- gene utilities --------------
def _resolve_gene_name(ad: sc.AnnData, query: Optional[str]) -> str:
    var = pd.Index(ad.var_names)
    lower_map = {v.lower(): v for v in var}
    if query:
        q = query.strip().lower()
        if q in lower_map:
            return lower_map[q]
        hits = [v for v in var if q in v.lower()]
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            print("\nMultiple genes matched your query:")
            for i, h in enumerate(hits[:50]):
                print(f"{i:>3}  {h}")
            sel = input("Pick index from matches (or type name): ").strip()
            try:
                i = int(sel); return hits[i]
            except Exception:
                return sel if sel in var else hits[0]
    print("\nNo gene specified. Example genes:\n", ", ".join(list(var[:20])), "...")
    g = input("Type gene symbol exactly as in var_names: ").strip()
    if g in ad.var_names:
        return g
    if g.lower() in lower_map:
        return lower_map[g.lower()]
    raise SystemExit(f"Gene '{g}' not found in ad.var_names.")

def _get_expr(ad: sc.AnnData, gene: str, source: str = "X", log1p: bool = False) -> np.ndarray:
    source = (source or "X").strip().lower()
    if source == "raw":
        if ad.raw is None:
            raise SystemExit("ad.raw is None, but EXPR_SOURCE='raw'.")
        x = ad.raw[:, gene].X
    elif source == "x":
        x = ad[:, gene].X
    else:
        if source not in ad.layers:
            raise SystemExit(f"Layer '{source}' not found in ad.layers.")
        x = ad[:, gene].layers[source]
    if hasattr(x, "toarray"):
        x = x.toarray().ravel()
    else:
        x = np.asarray(x).ravel()
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    if log1p:
        x = np.log1p(x)
    return x

def _compute_robust_limits(x: np.ndarray, vminp: Optional[float], vmaxp: Optional[float]) -> Tuple[float,float]:
    if vminp is None or vmaxp is None:
        return float(np.nanmin(x)), float(np.nanmax(x))
    lo = np.nanpercentile(x, float(vminp))
    hi = np.nanpercentile(x, float(vmaxp))
    if hi <= lo:
        hi = float(np.nanmax(x))
        lo = float(np.nanmin(x))
    return float(lo), float(hi)


# -------------- main --------------
def main():
    # (1) Load h5ad
    h5 = Path(SUBCLUSTER_H5AD) if SUBCLUSTER_H5AD else _pick_file("Select subcluster .h5ad", (".h5ad",".hdf5",".h5"))
    if not h5 or not h5.exists():
        raise SystemExit("No .h5ad selected.")
    ad = sc.read_h5ad(str(h5))
    df = ad.obs.copy()
    for k in ("spatial","X_spatial","spatial_global"):
        if k in ad.obsm and isinstance(ad.obsm[k], np.ndarray) and ad.obsm[k].shape[1] >= 2:
            df["X_obsm0"] = ad.obsm[k][:,0]
            df["Y_obsm1"] = ad.obsm[k][:,1]
            break

    # (2) columns
    region_col = next((k for k in REGION_COL_HINTS if k in df.columns), None)
    if region_col is None:
        cands = []
        for c in df.columns:
            s = df[c]
            if s.dtype.name in ("category","object") and 2 <= s.nunique() <= min(200, len(s)//10):
                if any(t in c.lower() for t in ["dataset","section","region","slide","fov","tile"]):
                    cands.append(c)
        region_col = cands[0] if cands else df.columns[0]
    label_col = next((k for k in LABEL_COL_HINTS if k in df.columns), None)

    # (3) choose regions & gene
    region_values = _choose_regions(df, region_col)
    gene = _resolve_gene_name(ad, GENE_NAME)
    expr_all = _get_expr(ad, gene, source=EXPR_SOURCE, log1p=LOG1P)

    # (4) DAPI
    dapi_path = Path(DAPI_FILE) if DAPI_FILE else _pick_file("Select the DAPI mosaic .tif for THIS section", (".tif",".tiff"))
    if not dapi_path or not dapi_path.exists():
        raise SystemExit("No DAPI mosaic selected.")
    img, (W_full, H_full), lvl = _read_pyramid_image(dapi_path)
    H_lvl, W_lvl = img.shape[0], img.shape[1]
    print(f"[DAPI] full={W_full}×{H_full}  level={lvl} → {W_lvl}×{H_lvl}")

    # (5) xy columns
    xcol, ycol = _find_xy_cols(df)
    x_src = df[xcol].to_numpy(dtype=float)
    y_src = df[ycol].to_numpy(dtype=float)

    # (6) transform handling
    if SKIP_TRANSFORM:
        x_full_all = x_src
        y_full_all = y_src
        coord_note = "pixel (skip transform)"
        transform_used = False
    else:
        unit = _infer_unit(x_src, y_src, W_full, H_full)

        print("\n(Optional) If you have a micron→mosaic pixel transform (e.g., experiment.json), pick it now.")
        tpath = _pick_file("Select transform file (Cancel to skip)", (".json",".txt",".csv",".vzg",".vzg2"))
        if tpath and tpath.is_dir():
            tpath = None

        if tpath:
            M = _load_transform_matrix(tpath)
            if FORCE_TRANSFORM_DIR in ("forward", "inverse"):
                if FORCE_TRANSFORM_DIR == "forward":
                    xf, yf = _apply_affine(M, x_src, y_src)
                    how, frac = "forward(force)", None
                else:
                    from numpy.linalg import inv
                    Mi = inv(M)
                    xf, yf = _apply_affine(Mi, x_src, y_src)
                    how, frac = "inverse(force)", None
                x_full_all, y_full_all = xf, yf
                coord_note = f"affine:{how}"
            else:
                how, score, xf, yf, frac = _auto_transform_direction(M, x_src, y_src, W_full, H_full)
                x_full_all, y_full_all = xf, yf
                coord_note = f"affine:{how}"
                print(f"[transform] chose {coord_note} (inside-frac ~ {frac:.2f})")
            transform_used = True
        else:
            # no transform picked → rely on units
            if unit == "micron":
                if UM_PER_PIXEL_FALLBACK is None:
                    raise SystemExit("Micron-looking coords but UM_PER_PIXEL_FALLBACK is None.")
                x_full_all = x_src / float(UM_PER_PIXEL_FALLBACK)
                y_full_all = y_src / float(UM_PER_PIXEL_FALLBACK)
                coord_note = f"micron→pixel / {UM_PER_PIXEL_FALLBACK} µm/px"
            else:
                x_full_all = x_src
                y_full_all = y_src
                coord_note = "pixel (no transform)"
            transform_used = False

    # (7) scale to chosen pyramid level
    xs_all, ys_all = _scale_to_level(x_full_all, y_full_all, W_full, H_full, W_lvl, H_lvl)
    if INVERT_Y:
        ys_all = H_lvl - ys_all

    # QC
    print(f"[qc] xs range: {np.nanmin(xs_all):.1f}..{np.nanmax(xs_all):.1f}  "
          f"ys range: {np.nanmin(ys_all):.1f}..{np.nanmax(ys_all):.1f}  "
          f"img level size: {W_lvl}×{H_lvl}")

    # (8) attach arrays to df
    df["_xs_"] = xs_all
    df["_ys_"] = ys_all
    df["_expr_"] = expr_all
    
    # ---- Global color scale (once for all selected regions) ----
    if GLOBAL_COLOR_SCALE:
        # mask rows that belong to any of the chosen regions
        mask = df[region_col].astype(str).isin([str(r) for r in region_values])
        expr_sel = df.loc[mask, "_expr_"].to_numpy(dtype=float)
    
        # allow locking to absolute numbers if user sets FIXED_VMIN/FIXED_VMAX
        if FIXED_VMIN is not None and FIXED_VMAX is not None:
            global_vmin, global_vmax = float(FIXED_VMIN), float(FIXED_VMAX)
        else:
            global_vmin, global_vmax = _compute_robust_limits(expr_sel, VMIN_PERC, VMAX_PERC)
            if not np.isfinite(global_vmin): global_vmin = float(np.nanmin(expr_sel))
            if not np.isfinite(global_vmax): global_vmax = float(np.nanmax(expr_sel))
            if global_vmax <= global_vmin: global_vmax = global_vmin + 1e-9
    
        print(f"[scale] using GLOBAL color scale v∈[{global_vmin:.3g},{global_vmax:.3g}]")
    else:
        global_vmin = global_vmax = None  # computed per-section later

    if PLOT_CATEGORIES_OVER and label_col is None:
        print("[note] PLOT_CATEGORIES_OVER=True but no label column found; disabling overlay.")
        overlay_labels = None
    else:
        overlay_labels = label_col

    # (9) output dir
    out_dir = Path(OUT_DIR); out_dir.mkdir(parents=True, exist_ok=True)

    # (10) per-region plotting
    for region_value in region_values:
        sub = df[df[region_col].astype(str) == str(region_value)].copy()
        if sub.empty:
            sub = df[df[region_col].astype(str).str.contains(str(region_value), case=False, na=False)].copy()
        if sub.empty:
            print(f"[skip] No rows for region '{region_value}' in column '{region_col}'")
            continue

        xs = sub["_xs_"].to_numpy(dtype=float)
        ys = sub["_ys_"].to_numpy(dtype=float)
        expr = sub["_expr_"].to_numpy(dtype=float)

        # Downsample
        if MAX_POINTS and xs.size > MAX_POINTS:
            idx = np.random.choice(xs.size, MAX_POINTS, replace=False)
            xs, ys, expr = xs[idx], ys[idx], expr[idx]
            if overlay_labels is not None:
                labs = sub[overlay_labels].astype(str).to_numpy()[idx]
            print(f"[downsample] plotting {len(xs):,} of {len(sub):,} cells in '{region_value}'")
        else:
            if overlay_labels is not None:
                labs = sub[overlay_labels].astype(str).to_numpy()

        # robust range
        # color scale
        if GLOBAL_COLOR_SCALE:
            vmin, vmax = global_vmin, global_vmax
        else:
            vmin, vmax = _compute_robust_limits(expr, VMIN_PERC, VMAX_PERC)
            if not np.isfinite(vmin): vmin = float(np.nanmin(expr))
            if not np.isfinite(vmax): vmax = float(np.nanmax(expr))
            if vmax <= vmin: vmax = vmin + 1e-9

        # inside fraction
        inside = (xs >= -5) & (xs <= W_lvl+5) & (ys >= -5) & (ys <= H_lvl+5)
        print(f"[sanity] '{region_value}': points inside image: {inside.mean():.2%}")

        # plot
        safe_region = re.sub(r'[^A-Za-z0-9_.-]+','_',str(region_value))
        out_file = out_dir / f"overlay_gene-{re.sub(r'[^A-Za-z0-9_.-]+','_',gene)}__{safe_region}.svg"

        plt.figure(figsize=(min(12, max(6, W_lvl/1000)), min(12, max(6, H_lvl/1000))), dpi=DPI)
        plt.imshow(img, cmap="gray", origin="upper")
        sc1 = plt.scatter(xs, ys, s=DOT_SIZE, c=np.clip(expr, vmin, vmax),
                          vmin=vmin, vmax=vmax, cmap=CMAP, alpha=ALPHA, linewidths=0, edgecolors="none")
        plt.xlim([0, W_lvl]); plt.ylim([H_lvl, 0]); plt.xticks([]); plt.yticks([])

        ttl = (f"{Path(dapi_path).name} | region='{Path(str(region_value)).name}' | n={len(expr):,} | "
               f"{coord_note} | gene={gene} [{EXPR_SOURCE}{' +log1p' if LOG1P else ''}] "
               f"| v∈[{vmin:.3g},{vmax:.3g}]")
        plt.title(ttl, fontsize=9)

        cbar = plt.colorbar(sc1, fraction=0.046, pad=0.04)
        cbar.set_label(f"{gene} expression", rotation=90, fontsize=8)

        # optional faint cluster overlay
        if PLOT_CATEGORIES_OVER and overlay_labels is not None:
            cats = sorted(pd.unique(labs))
            if len(cats) <= MAX_LEGEND:
                palette = _choose_palette(len(cats))
                lut = {c: palette[i] for i, c in enumerate(cats)}
                cols = np.array([lut[c] for c in labs])
                plt.scatter(xs, ys, s=DOT_SIZE*0.6, c=cols, alpha=0.25, linewidths=0, edgecolors="none")
                if LEGEND:
                    from matplotlib.lines import Line2D
                    handles = [Line2D([0],[0], marker='o', linestyle='', markersize=6,
                                      markerfacecolor=lut[c], label=c, alpha=0.6) for c in cats]
                    plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5),
                               title=overlay_labels, frameon=False, fontsize=7)

        plt.tight_layout()
        plt.savefig(out_file, dpi=DPI)
        plt.close()
        print(f"[SAVED] {out_file}")

if __name__ == "__main__":
    main()
