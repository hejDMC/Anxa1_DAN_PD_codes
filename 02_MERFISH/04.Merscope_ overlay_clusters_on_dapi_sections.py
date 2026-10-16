from __future__ import annotations
import os, re, json, math
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import scanpy as sc
import tifffile

# ============ USER TWEAKS ============
SUBCLUSTER_H5AD: Optional[str] = r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad"
REGION_COL_HINTS = ["Datasets", "section_id", "Region", "region", "slide_id"]
LABEL_COL_HINTS  = ["leiden_sub_Th_Slc6a3", "leiden_sub", "DA_subclusters_Th_Slc6a3", "leiden"]

START_DAPI_DIR = r"Y:\Vasska\Merscope_data\202507171430_Mitopark16w2ko_VMSC13902\region_ko1c\images"
DAPI_FILE      = r"Y:\Vasska\Merscope_data\202507171430_Mitopark16w2ko_VMSC13902\region_ko1c\images\mosaic_DAPI_z3.tif"

PYRAMID_LEVEL = -1      # 0=full, -1=lowest
DOT_SIZE = 10
ALPHA = 0.85
LEGEND = True
MAX_LEGEND = 25
DPI = 180
OUT_DIR = r"F:/vizgen_clustering_out/overlays_one_section/each_cluster"

UM_PER_PIXEL_FALLBACK: Optional[float] = 0.108
INVERT_Y = False
MAX_POINTS = 300_000
# =====================================

# ---- tiny GUI pickers ----
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

# ---- helpers ----
def _find_xy_cols(df: pd.DataFrame) -> Tuple[str,str]:
    cands = [("global_x","global_y"),("x_global","y_global"),("x","y"),
             ("center_x","center_y"),("x_centroid","y_centroid")]
    lower = {c.lower(): c for c in df.columns}
    for a,b in cands:
        if a in lower and b in lower: return lower[a], lower[b]
    for a,b in [("X_obsm0","Y_obsm1")]:
        if a in df and b in df: return a,b
    raise KeyError("Could not find XY columns (looked for global_x/global_y, x_global/y_global, center_x/center_y, x/y).")

def _choose_region(df: pd.DataFrame, region_col: str) -> str:
    vals = sorted(df[region_col].astype(str).unique().tolist())
    print("\nAvailable regions:")
    for i, v in enumerate(vals):
        print(f"{i:>3}  {v}")
    sel = input("\nType index of the region you want (e.g., 0): ").strip()
    try:
        i = int(sel)
        return vals[i]
    except Exception:
        hit = [v for v in vals if sel and sel.lower() in v.lower()]
        return hit[0] if hit else vals[0]

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
        keys = ["micron_to_mosaic_pixel_transform","micron_to_mosaic_pixel_matrix","microns_to_mosaic_pixels",
                "mosaic_to_micron","mosaic_pixel_to_micron_transform"]
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
    flat = np.loadtxt(path, dtype=float).ravel()
    if flat.size == 6:  M = flat.reshape(2,3)
    elif flat.size == 9: M = flat.reshape(3,3)
    else: raise ValueError("Transform text should have 6 or 9 numbers.")
    return _augment_to_3x3(M)

def _apply_affine(M: np.ndarray, x_um: np.ndarray, y_um: np.ndarray) -> Tuple[np.ndarray,np.ndarray]:
    P = np.vstack([x_um, y_um, np.ones_like(x_um, float)])
    Q = (M @ P).T
    return Q[:,0], Q[:,1]

def _auto_transform_direction(M: np.ndarray, x_um: np.ndarray, y_um: np.ndarray, W: int, H: int):
    from numpy.linalg import inv
    xf, yf = _apply_affine(M, x_um, y_um)
    fi = (xf >= -0.05*W) & (xf <= 1.05*W) & (yf >= -0.05*H) & (yf <= 1.05*H)
    frac_fwd = fi.mean()
    try:
        Mi = inv(M)
        xi, yi = _apply_affine(Mi, x_um, y_um)
        fi2 = (xi >= -0.05*W) & (xi <= 1.05*W) & (yi >= -0.05*H) & (yi <= 1.05*H)
        frac_inv = fi2.mean()
    except Exception:
        frac_inv = -1.0
    if frac_inv > frac_fwd:
        return "inverse", xi, yi, frac_inv
    else:
        return "forward", xf, yf, frac_fwd

def _infer_unit(x: np.ndarray, y: np.ndarray, W: int, H: int) -> str:
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if x.size == 0: return "pixel"
    xq, yq = np.quantile(np.abs(x), 0.95), np.quantile(np.abs(y), 0.95)
    return "pixel" if (xq <= W*1.5 and yq <= H*1.5) else "micron"

# ---- main ----
def main():
    # 1) subcluster object
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

    # 2) columns
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
    if label_col is None:
        raise SystemExit("Could not find a subcluster label column (looked for leiden_sub_Th_Slc6a3, leiden_sub, DA_subclusters_Th_Slc6a3, leiden).")

    # 3) pick region
    region_value = _choose_region(df, region_col)
    sub = df[df[region_col].astype(str) == str(region_value)].copy()
    if sub.empty:
        sub = df[df[region_col].astype(str).str.contains(str(region_value), case=False, na=False)].copy()
    if sub.empty:
        raise SystemExit(f"No rows for region '{region_value}' in column '{region_col}'")
    print(f"\n[region] {region_col}='{region_value}'  → cells: {len(sub):,}")
    
    # 4) Choose the cluster to map
    print("\nAvailable clusters:")
    unique_clusters = sorted(df[label_col].unique())
    for i, cluster in enumerate(unique_clusters):
        print(f"{i:>3}  {cluster}")
    cluster_idx = input("\nType the index of the cluster you want to map (e.g., 0): ").strip()
    
    try:
        cluster_idx = int(cluster_idx)
        selected_cluster = unique_clusters[cluster_idx]
    except ValueError:
        selected_cluster = cluster_idx
        if selected_cluster not in unique_clusters:
            print(f"Invalid cluster, defaulting to {unique_clusters[0]}")
            selected_cluster = unique_clusters[0]
    
    # Filter the data to include only the selected cluster
    sub = sub[sub[label_col] == selected_cluster]
    print(f"[selected cluster] {selected_cluster} → cells: {len(sub):,}")

    # 4) DAPI
    # 4) DAPI
    dapi_path = Path(DAPI_FILE) if DAPI_FILE else _pick_file("Select the DAPI mosaic .tif for THIS section", (".tif",".tiff"))
    if not dapi_path or not dapi_path.exists():
        print("[Debug] No DAPI mosaic selected or file doesn't exist. Exiting.")
        raise SystemExit("No DAPI mosaic selected.")
    else:
        img, (W_full, H_full), lvl = _read_pyramid_image(dapi_path)
        H_lvl, W_lvl = img.shape[0], img.shape[1]
        print(f"[DAPI] full={W_full}×{H_full}  level={lvl} → {W_lvl}×{H_lvl}")

    # 5) XY
    xcol, ycol = _find_xy_cols(sub)
    x = sub[xcol].astype(float).to_numpy()
    y = sub[ycol].astype(float).to_numpy()

    # 6) optional transform
    print("\n(Optional) If you have a micron→mosaic pixel transform (e.g., experiment.json), pick it now.")
    tpath = _pick_file("Select transform file (Cancel to skip)", (".json",".txt",".csv",".vzg",".vzg2"))
    if tpath and tpath.is_dir():
        tpath = None

    if tpath:
        M = _load_transform_matrix(tpath)
        how, x_full, y_full, frac = _auto_transform_direction(M, x, y, W_full, H_full)
        print(f"[transform] chose {how} mapping (inside-frac ~ {frac:.2f})")
        coord_note = "affine"
    else:
        unit = _infer_unit(x, y, W_full, H_full)
        if unit == "micron":
            if UM_PER_PIXEL_FALLBACK is None:
                raise SystemExit("Coords look like microns but UM_PER_PIXEL_FALLBACK is None. Set a value (e.g., 0.108).")
            x_full, y_full = x / float(UM_PER_PIXEL_FALLBACK), y / float(UM_PER_PIXEL_FALLBACK)
            coord_note = f"micron→pixel / {UM_PER_PIXEL_FALLBACK} µm/px"
        else:
            x_full, y_full = x, y
            coord_note = "pixel (no transform)"

    # 7) scale to level
    xs, ys = _scale_to_level(x_full, y_full, W_full, H_full, W_lvl, H_lvl)
    if INVERT_Y:
        ys = H_lvl - ys

    # 8) color by subcluster
    labels = sub[label_col].astype(str).to_numpy()
    if MAX_POINTS and xs.size > MAX_POINTS:
        idx = np.random.choice(xs.size, MAX_POINTS, replace=False)
        xs, ys, labels = xs[idx], ys[idx], labels[idx]
        print(f"[downsample] plotting {len(xs):,} of {len(sub):,} cells")

    cats = sorted(pd.unique(labels))
    palette = _choose_palette(len(cats))
    lut = {c: palette[i] for i, c in enumerate(cats)}
    cols = np.array([lut[c] for c in labels])

    # sanity
    inside = (xs >= -5) & (xs <= W_lvl+5) & (ys >= -5) & (ys <= H_lvl+5)
    print(f"[sanity] points inside image: {inside.mean():.2%}")

    # 9) plot & save
    out_dir = Path(OUT_DIR); out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"overlay__{re.sub(r'[^A-Za-z0-9_.-]+','_',str(region_value))}.svg"

    plt.figure(figsize=(min(12, max(6, W_lvl/1000)), min(12, max(6, H_lvl/1000))), dpi=DPI)
    plt.imshow(img, cmap="gray", origin="upper")
    plt.scatter(xs, ys, s=DOT_SIZE, c=cols, alpha=ALPHA, linewidths=0, edgecolors="none")
    plt.xlim([0, W_lvl]); plt.ylim([H_lvl, 0]); plt.xticks([]); plt.yticks([])
    ttl = f"{Path(dapi_path).name}  |  region='{Path(str(region_value)).name}'  |  n={len(labels):,}  |  {coord_note}"
    plt.title(ttl, fontsize=9)

    if LEGEND and len(cats) <= MAX_LEGEND:
        from matplotlib.lines import Line2D
        handles = [Line2D([0],[0], marker='o', linestyle='', markersize=6,
                          markerfacecolor=lut[c], label=c) for c in cats]
        plt.legend(handles=handles, loc='center left', bbox_to_anchor=(1.02, 0.5),
                   title=label_col, frameon=False, fontsize=7)

    plt.tight_layout()
    plt.savefig(out_file, dpi=DPI)
    plt.close()
    print(f"\n[SAVED] {out_file}\n")

if __name__ == "__main__":
    main()
