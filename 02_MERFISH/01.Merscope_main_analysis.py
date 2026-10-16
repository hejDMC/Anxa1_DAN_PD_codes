# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 12:32:38 2025

@author: Vasiliki Skara

"""

#%% Load MERSCOPE data and inspect metadata
from __future__ import annotations
import numpy as np, pandas as pd, scanpy as sc

# If you don't already have 'adata' in memory:
adata = sc.read_h5ad(r"F:/14-10-25_10-07_All_animals_exp_Cells.hdf5")

def summarize_anndata(adata, n=20):
    print("=== AnnData ===")
    print(adata)
    print("\n=== .obs (first 50 columns) ===")
    print(list(adata.obs.columns)[:50])
    print("\n=== .var shape & sample gene names ===")
    print(adata.var.shape, list(adata.var_names[:n]))
    print("\n=== slots ===")
    print("obsm:", list(adata.obsm.keys()))
    print("layers:", list(adata.layers.keys()))
    print("uns:", list(adata.uns.keys()))

    # Heuristic: find obs columns likely to be batch indicators
    candidates = []
    for c in adata.obs.columns:
        cname = c.lower()
        if any(k in cname for k in ["dataset","section","slide","animal","batch","group"]):
            candidates.append(c)
    # Also add low-cardinality categorical-like columns
    for c in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[c]) or pd.api.types.is_categorical_dtype(adata.obs[c]):
            nunq = adata.obs[c].nunique(dropna=True)
            if 1 < nunq <= min(50, max(2, adata.n_obs // 100)):
                if c not in candidates:
                    candidates.append(c)

    print("\n=== candidate batch columns ===")
    print(candidates or "(none found)")
    for c in candidates:
        vc = adata.obs[c].astype(str).value_counts(dropna=False)
        print(f"\n-- {c} (n_unique={vc.shape[0]}) --")
        print(vc.head(15))

summarize_anndata(adata)

# If 'Datasets' exists, promote it to a batch key used by the pipeline:
if 'Datasets' in adata.obs:
    adata.obs['section_id'] = adata.obs['Datasets'].astype(str)
    print("\nCreated adata.obs['section_id'] from 'Datasets'.")

# If you'd rather use your manual grouping:
# if 'Custom Cells groups' in adata.obs:
#     adata.obs['section_id'] = adata.obs['Custom Cells groups'].astype(str)
#     print("\nCreated adata.obs['section_id'] from 'Custom Cells groups'.")

# Quick visual check (will make sense after running neighbors/UMAP):
if 'section_id' in adata.obs:
    try:
        sc.pl.umap(adata, color=['section_id'], wspace=0.4)
    except Exception as e:
        print("UMAP not present yet — run neighbors/UMAP first. Error:", e)
        

#%% Quality control, normalization, and whole-dataset clustering
import os
from glob import glob
from pathlib import Path
from typing import Optional, Iterable, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse as sp

# ------------------ CONFIG (edit these) ------------------ #
# If None, will auto-glob the first '*.hdf5' not ending with '_clustered.hdf5'
ANNDATA_PATH: Optional[str] =  r"F:/14-10-25_10-07_All_animals_exp_Cells.hdf5"  

if 'section_id' not in adata.obs and 'Datasets' in adata.obs:
    adata.obs['section_id'] = adata.obs['Datasets'].astype(str)

# Clustering params (optimized defaults for targeted panels)
N_NEIGHBORS = 12
N_PCS = 30
LEIDEN_RESOLUTION = 1
UMAP_MIN_DIST = 0.3
NEIGHBOR_METRIC = "cosine"   # cosine works well for composition-like profiles
RANDOM_STATE = 42

# Batch integration (BBKNN). Either force with USE_BBKNN=True
# or let AUTO_BBKNN decide based on diagnostics.
USE_BBKNN = False
AUTO_BBKNN = True
BATCH_KEY: Optional[str] = None     # e.g., 'section_id', 'slide_id'
BATCH_KEY_PRIORITY = ['section_id','Datasets','slide_id','animal_id']
BBKNN_NEIGHBORS_WITHIN_BATCH: Optional[int] = None  # if None → max(3, N_NEIGHBORS//3)
POST_BBKNN_LEIDEN_RESOLUTION = LEIDEN_RESOLUTION    # e.g., set to 0.9 if you want less splitting
MIXING_SCORE_MIN = 0.30   # if neighbor-mixing < this → integrate
CLUSTER_PURITY_MAX = 0.70 # if many clusters have purity > this → integrate

# QC / filtering (set to None to skip a criterion)
APPLY_QC_FILTERS = True
MIN_COUNTS_PER_CELL: Optional[int] = 25   # Vizgen example: total expression > 25
MIN_GENES_PER_CELL: Optional[int] = None  # often not used for panels
MAX_PCT_MT: Optional[float] = None        # e.g., 20.0

# Remove Blank probes? (recommended)
REMOVE_BLANKS = True
BLANK_KEYWORD = 'Blank'  # probes containing this substring in var_names will be dropped

# HVGs (OFF by default for panels; the template originally used all genes)
USE_HVG = False
HVG_N_TOP_GENES = 2000
HVG_FLAVOR = 'seurat_v3'

# Scaling & optional regression
DO_SCALE = True
SCALE_MAX_VALUE = 10
REGRESS_OUT_TOTAL_COUNTS = False   # set True if total_counts drives the embedding

# Diagnostics & utilities
OUT_DIR = Path("F:/vizgen_clustering_out")
SAVE_PLOTS = True
PLOT_DPI = 150
RUN_PARAM_SWEEP = False  # try a tiny grid of (neighbors, pcs, min_dist, resolution)

# --------------------------------------------------------- #

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=PLOT_DPI)
np.random.seed(RANDOM_STATE)


# ---------- Helpers ---------- #

def _pick_input_path() -> Path:
    if ANNDATA_PATH:
        return Path(ANNDATA_PATH)
    candidates = [x for x in glob('*.hdf5') if not x.endswith('_clustered.hdf5')]
    if not candidates:
        candidates = [x for x in glob('*.h5ad') if not x.endswith('_clustered.h5ad')]
    if not candidates:
        raise FileNotFoundError("No input .hdf5/.h5ad found. Set ANNDATA_PATH or place the file next to this script.")
    print(f"Auto-selected input: {candidates[0]}")
    return Path(candidates[0])


def _set_mt_flag(adata: sc.AnnData) -> None:
    """Create adata.var['mt'] boolean for mitochondrial genes (handles 'MT-'/'mt-')."""
    gene_names = np.asarray(adata.var_names).astype(str)
    mt_mask = np.char.startswith(np.char.upper(gene_names), 'MT-')
    adata.var['mt'] = pd.Series(mt_mask, index=adata.var_names)


def _qc_metrics(adata: sc.AnnData) -> None:
    """Compute QC metrics safely for small targeted panels (no percent_top)."""
    _set_mt_flag(adata)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)


def _apply_filters(adata: sc.AnnData) -> sc.AnnData:
    if not APPLY_QC_FILTERS:
        return adata
    mask = np.ones(adata.n_obs, dtype=bool)
    if MIN_COUNTS_PER_CELL is not None:
        mask &= adata.obs['total_counts'].to_numpy() >= MIN_COUNTS_PER_CELL
    if MIN_GENES_PER_CELL is not None and 'n_genes_by_counts' in adata.obs:
        mask &= adata.obs['n_genes_by_counts'].to_numpy() >= MIN_GENES_PER_CELL
    if MAX_PCT_MT is not None and 'pct_counts_mt' in adata.obs:
        mask &= adata.obs['pct_counts_mt'].to_numpy() <= MAX_PCT_MT
    before = adata.n_obs
    adata = adata[mask].copy()
    print(f"Filtered cells: {before} → {adata.n_obs}")
    return adata


def _remove_blanks(adata: sc.AnnData) -> sc.AnnData:
    if not REMOVE_BLANKS:
        return adata
    var = pd.Index(adata.var_names.astype(str))
    keep = ~var.str.contains(BLANK_KEYWORD, case=False, regex=False)
    before = adata.n_vars
    adata = adata[:, keep].copy()
    print(f"Removed blanks: genes {before} → {adata.n_vars}")
    return adata


def _maybe_hvg(adata: sc.AnnData) -> None:
    if not USE_HVG:
        return
    n_top = min(HVG_N_TOP_GENES, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor=HVG_FLAVOR, subset=False, inplace=True)
    hvgs = adata.var['highly_variable'] if 'highly_variable' in adata.var else None
    if hvgs is not None:
        print(f"Selected HVGs: {int(hvgs.sum())}")


def _cluster_pipeline(adata: sc.AnnData) -> None:
    # Normalize & log1p
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Optional: regress out library size (can sharpen manifolds when counts vary)
    if REGRESS_OUT_TOTAL_COUNTS and 'total_counts' in adata.obs:
        sc.pp.regress_out(adata, keys=['total_counts'])

    # Optional HVGs (template originally did not do HVG selection)
    _maybe_hvg(adata)

    # Scale before PCA; helps targeted panels with uneven dynamic ranges
    if DO_SCALE:
        sc.pp.scale(adata, max_value=SCALE_MAX_VALUE)

    # PCA (use HVGs only if enabled)
    sc.tl.pca(adata, svd_solver='arpack', use_highly_variable=USE_HVG, n_comps=N_PCS, random_state=RANDOM_STATE)

    # Graph: start with standard neighbors; optionally replace with BBKNN based on config/diagnostics
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, metric=NEIGHBOR_METRIC)
    _maybe_run_bbknn(adata)

    # UMAP + Leiden
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, random_state=RANDOM_STATE)
    res = POST_BBKNN_LEIDEN_RESOLUTION if adata.uns.get('_used_bbknn', False) else LEIDEN_RESOLUTION
    sc.tl.leiden(adata, resolution=res)


def _export_vizgen_csv(adata: sc.AnnData, out_csv: Path) -> None:
    # Try to find a stable cell id column; fall back to index
    for k in ['cell_id', 'Cell_ID', 'cellID', 'CellId', 'cellId', 'cell_index']:
        if k in adata.obs:
            cell_ids = adata.obs[k].astype(str).values
            break
    else:
        cell_ids = adata.obs_names.astype(str).values

    df = pd.DataFrame({'cell_id': cell_ids, 'cluster': adata.obs['leiden'].astype(str).values})
    if 'X_umap' in adata.obsm:
        umap = pd.DataFrame(adata.obsm['X_umap'], columns=['umap_1','umap_2'])
        df = pd.concat([df.reset_index(drop=True), umap.reset_index(drop=True)], axis=1)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)
    print(f"Wrote Vizualizer CSV → {out_csv}")


def _save_diagnostics(adata: sc.AnnData, out_dir: Path) -> None:
    if not SAVE_PLOTS:
        return
    import matplotlib.pyplot as plt
    out_dir.mkdir(parents=True, exist_ok=True)

    sc.pl.umap(adata, color=['leiden'], save=None, show=False)
    plt.gcf().savefig(out_dir / 'umap_leiden.png', dpi=PLOT_DPI, bbox_inches='tight')

    cols = [c for c in ['total_counts','n_genes_by_counts','pct_counts_mt'] if c in adata.obs]
    if cols:
        sc.pl.umap(adata, color=cols, save=None, show=False)
        plt.gcf().savefig(out_dir / 'umap_qc.png', dpi=PLOT_DPI, bbox_inches='tight')

    batch_cols = [c for c in ['slide_id','section_id','animal_id','batch','Datasets'] if c in adata.obs]
    if batch_cols:
        sc.pl.umap(adata, color=batch_cols, save=None, show=False)
        plt.gcf().savefig(out_dir / 'umap_batch.png', dpi=PLOT_DPI, bbox_inches='tight')

    if 'pca' in adata.uns and 'variance_ratio' in adata.uns['pca']:
        sc.pl.pca_variance_ratio(adata, log=True, show=False)
        plt.gcf().savefig(out_dir / 'pca_variance_ratio.png', dpi=PLOT_DPI, bbox_inches='tight')

    bk = _choose_batch_key(adata)
    try:
        with open(out_dir / 'batch_diagnostics.txt', 'w') as fh:
            fh.write(f"used_bbknn: {adata.uns.get('_used_bbknn', False)}\n")
            if adata.uns.get('_used_bbknn', False):
                fh.write(f"bbknn_batch_key: {adata.uns.get('_bbknn_batch_key', bk)}\n")
            if bk is not None:
                mix = _neighbor_mixing_score(adata, bk)
                fh.write(f"neighbor_mixing({bk}) = {mix}\n")
                if 'leiden' in adata.obs:
                    purity = _cluster_batch_purity(adata, 'leiden', bk)
                    fh.write("cluster_purity_top10:\n")
                    for idx, val in purity.sort_values(ascending=False).head(10).items():
                        fh.write(f"  {idx}: {val:.3f}\n")
    except Exception as e:
        print(f"Batch diagnostics write failed: {e}")



def _silhouette_param_sweep(adata: sc.AnnData, grid: Iterable[Tuple[int,int,float,float]]) -> Tuple[Tuple[int,int,float,float], float]:
    from sklearn.metrics import silhouette_score
    best, best_score = None, -1
    for nn, pcs, res, mind in grid:
        sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=pcs, metric=NEIGHBOR_METRIC)
        sc.tl.umap(adata, min_dist=mind, random_state=RANDOM_STATE)
        sc.tl.leiden(adata, resolution=res, key_added='leiden_tmp')
        emb = adata.obsm['X_umap']
        labels = adata.obs['leiden_tmp'].astype('category').cat.codes.to_numpy()
        if emb.shape[0] > 10000:
            # sample for speed
            idx = np.random.choice(emb.shape[0], 10000, replace=False)
            s = silhouette_score(emb[idx], labels[idx])
        else:
            s = silhouette_score(emb, labels)
        print(f"grid nn={nn}, pcs={pcs}, res={res}, mind={mind} → silhouette={s:.4f}")
        if s > best_score:
            best, best_score = (nn, pcs, res, mind), s
    return best, best_score


def hdf5_tree(path: str) -> None:
    """Print the HDF5/AnnData group/dataset tree (for exploration)."""
    import h5py
    def walk(g, p=""):
        for k, v in g.items():
            q = f"{p}/{k}"
            typ = 'DATASET' if isinstance(v, h5py.Dataset) else 'GROUP'
            print(f"[{typ}] {q}")
            if isinstance(v, h5py.Group):
                walk(v, q)
    with h5py.File(path, 'r') as f:
        walk(f)


def _choose_batch_key(adata: sc.AnnData) -> Optional[str]:
    """Pick a batch key if not explicitly set, based on presence in .obs.
    Priority comes from BATCH_KEY_PRIORITY.
    """
    if BATCH_KEY is not None and BATCH_KEY in adata.obs.columns:
        return BATCH_KEY
    for k in BATCH_KEY_PRIORITY:
        if k in adata.obs.columns:
            return k
    return None


def _neighbor_mixing_score(adata: sc.AnnData, batch: str) -> float:
    """1 - mean(same-batch fraction among neighbors). Higher is better mixing.
    Requires a neighbors graph in adata.obsp.
    """
    key = adata.uns.get('neighbors', {}).get('connectivities_key', 'connectivities')
    if key not in adata.obsp:
        return np.nan
    A = adata.obsp[key].tocsr()
    b = adata.obs[batch].astype('category').cat.codes.to_numpy()
    same = []
    for i in range(A.shape[0]):
        nbrs = A.indices[A.indptr[i]:A.indptr[i+1]]
        if len(nbrs) == 0:
            continue
        same.append((b[nbrs] == b[i]).mean())
    return 1.0 - float(np.mean(same)) if same else np.nan


def _cluster_batch_purity(adata: sc.AnnData, group: str, batch: str) -> pd.Series:
    """Max fraction of any batch within each cluster (higher → more batchy)."""
    ct = pd.crosstab(adata.obs[group], adata.obs[batch])
    prop = ct.div(ct.sum(1), 0)
    return prop.max(1)


def _maybe_run_bbknn(adata: sc.AnnData) -> None:
    """Decide whether to run BBKNN based on config and diagnostics.
    Sets adata.uns['_used_bbknn']=True if BBKNN was applied.
    """
    adata.uns['_used_bbknn'] = False

    # If user forced it, just run
    batch_key = _choose_batch_key(adata)
    if USE_BBKNN and batch_key is not None:
        try:
            import bbknn
            nwb = BBKNN_NEIGHBORS_WITHIN_BATCH or max(3, N_NEIGHBORS // 3)
            bbknn.bbknn(adata, batch_key=batch_key, neighbors_within_batch=nwb, n_pcs=N_PCS)
            adata.uns['_used_bbknn'] = True
            adata.uns['_bbknn_batch_key'] = batch_key
            return
        except Exception as e:
            print(f"BBKNN skipped (import/run error): {e}")

    if not AUTO_BBKNN or batch_key is None:
        return

    # Diagnostics on current graph
    mix = _neighbor_mixing_score(adata, batch_key)
    purity = _cluster_batch_purity(adata, group='leiden' if 'leiden' in adata.obs else None or 'tmp_cluster', batch=batch_key) if 'leiden' in adata.obs else pd.Series(dtype=float)

    should_integrate = False
    if not np.isnan(mix) and mix < MIXING_SCORE_MIN:
        should_integrate = True
    if not purity.empty and (purity > CLUSTER_PURITY_MAX).sum() >= max(1, len(purity)//5):
        should_integrate = True

    if should_integrate:
        try:
            import bbknn
            nwb = BBKNN_NEIGHBORS_WITHIN_BATCH or max(3, N_NEIGHBORS // 3)
            print(f"Running BBKNN (batch_key='{batch_key}', neighbors_within_batch={nwb}) due to mixing={mix:.3f}.")
            bbknn.bbknn(adata, batch_key=batch_key, neighbors_within_batch=nwb, n_pcs=N_PCS)
            adata.uns['_used_bbknn'] = True
            adata.uns['_bbknn_batch_key'] = batch_key
        except Exception as e:
            print(f"BBKNN auto-run failed: {e}")


# ---------- Main ---------- #

def run_pipeline(return_adata: bool = False) -> sc.AnnData:
    """Run the full pipeline and return the clustered AnnData.
    If called from Spyder, this also sets a global `adata` for convenience.
    """
    inp = _pick_input_path()
    print(f"Reading: {inp}")
    ad = sc.read_h5ad(str(inp))
    print(f"Loaded AnnData: {ad.n_obs} cells × {ad.n_vars} genes")

    # Remove blanks first (keeps memory down)
    ad = _remove_blanks(ad)

    # QC metrics + optional filters
    _qc_metrics(ad)
    ad = _apply_filters(ad)

    # Core clustering pipeline
    _cluster_pipeline(ad)

    # Optional: tiny parameter sweep (read-only; does not overwrite main results)
    if RUN_PARAM_SWEEP:
        grid = [(nn, pcs, res, md) for nn in (8,12,15) for pcs in (20,30) for res in (0.6,0.8,1.0) for md in (0.2,0.3,0.4)]
        best, best_score = _silhouette_param_sweep(ad, grid)
        print(f"Best grid = {best} with silhouette = {best_score:.4f}")

    # Save outputs
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    base = inp.with_suffix("")
    clustered = Path(str(base) + "_clustered.h5ad")
    ad.write_h5ad(clustered)
    print(f"Wrote clustered AnnData (.h5ad) → {clustered}")

    # Also save with .hdf5 extension for Vizgen/Vizualizer compatibility
    clustered_hdf5 = Path(str(base) + "_clustered.hdf5")
    ad.write_h5ad(clustered_hdf5)
    print(f"Wrote clustered AnnData (.hdf5) → {clustered_hdf5}")

    # Optional: save a compact per-cell HDF5 with just cluster annotations
    try:
        obs_cols = ['leiden'] + [c for c in ['section_id','slide_id','animal_id'] if c in ad.obs]
        obs_df = ad.obs[obs_cols].copy()
        obs_df.to_hdf(OUT_DIR / 'cluster_annotations.h5', key='clusters', mode='w')
        print(f"Wrote per-cell cluster annotations → {OUT_DIR / 'cluster_annotations.h5'}")
    except Exception as e:
        print(f"Saving cluster annotations HDF5 skipped: {e}")

    _save_diagnostics(ad, OUT_DIR)
    _export_vizgen_csv(ad, OUT_DIR / 'clusters_for_vizualizer.csv')

    # Expose to interactive session and/or return
    globals()['adata'] = ad
    return ad if return_adata else ad


def main():
    # Preserve behavior when running as a script
    run_pipeline(return_adata=False)


if __name__ == "__main__":
    main()

#%% Optional: compare Leiden clustering resolutions

import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path

def quick_leiden_resolution_sweep(
    adata,
    resolutions=(0.8, 0.9, 1.0, 1.2),
    out_dir: Path = OUT_DIR,
    point_size: float = 3.0,
    dpi: int = PLOT_DPI
):
    """
    Re-run Leiden at multiple resolutions on an already-neighbored/embedded AnnData.
    - Adds columns like 'leiden_r0.8' to .obs
    - Saves UMAP PNGs colored by each resolution
    - Prints cluster counts per resolution
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # Make sure the graph/UMAP exist (your pipeline already computes these)
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, metric=NEIGHBOR_METRIC)
    if 'X_umap' not in adata.obsm:
        sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, random_state=RANDOM_STATE)

    print("\n[resolution sweep]")
    for res in resolutions:
        key = f"leiden_r{res:g}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        nclu = adata.obs[key].nunique()
        print(f"  res={res:g} → {nclu} clusters (obs['{key}'])")

        sc.pl.umap(
            adata,
            color=[key],
            legend_loc="right margin",
            size=point_size,
            show=False
        )
        plt.gcf().savefig(out_dir / f"umap_{key}.png", dpi=dpi, bbox_inches="tight")
        plt.close('all')

    # Optional: also save a small summary CSV
    summary = (
        pd.Series({f"leiden_r{r:g}": adata.obs[f"leiden_r{r:g}"].nunique() for r in resolutions})
          .rename("n_clusters")
          .rename_axis("label")
          .to_frame()
    )
    summary.to_csv(out_dir / "leiden_resolution_summary.csv")
    print(f"Saved UMAPs and summary → {out_dir}")

# --- Usage ---
# ad = run_pipeline(return_adata=True)   # if you want to re-run from file
# OR just use your in-memory `adata` from the pipeline:
quick_leiden_resolution_sweep(adata, resolutions=(0.8, 0.9, 1.0, 1.2))

    


#%% Export whole-dataset UMAP and QC plots as PDFs

import pathlib, matplotlib.pyplot as plt, scanpy as sc
out = pathlib.Path(r"F:/vizgen_clustering_out"); out.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = str(out)                  # where Scanpy writes files when using `save=...`
sc.settings.file_format_figs = "pdf"           # default extension for `save=...`

# Make points rasterize inside PDFs (no need to pass rasterized=)
sc.settings.set_figure_params(dpi=150, dpi_save=1200, figsize=(10, 10), vector_friendly=True)
sc.settings.file_format_figs = "pdf"

ax = sc.pl.umap(
    adata,
    color='leiden',
    show=False,
    size=3,        # try 2–4 depending on density
    alpha=1.0      # avoid weird transparency in some PDF viewers
)
plt.gcf().savefig(out / "umap_leiden_rasterized.pdf", bbox_inches="tight", dpi=1200)
plt.close('all')


# QC overlays (only if present)
qc_cols = [c for c in ['total_counts','n_genes_by_counts','pct_counts_mt'] if c in adata.obs]
if qc_cols:
    sc.pl.umap(adata, color=qc_cols, show=False, save="_qc.pdf", s=3)

# Batch overlays (Datasets/section_id etc., if present)
batch_cols = [c for c in ['section_id','Datasets','slide_id','animal_id','batch'] if c in adata.obs]
if batch_cols:
    sc.pl.umap(adata, color=batch_cols, show=False, save="_batch.pdf", s=3)

# PCA variance ratio (if available from your run)
if 'pca' in adata.uns and 'variance_ratio' in adata.uns['pca']:
    sc.pl.pca_variance_ratio(adata, log=True, show=False, save="_pca_variance_ratio.pdf")

#%% Inspect animal distribution across the whole-dataset UMAP
# ==== UMAP colored by animal mapping (figure only; no file modifications) ====
import matplotlib.pyplot as plt
from pathlib import Path, PureWindowsPath
import pandas as pd
import numpy as np

# 1) Your mapping (as given)
dataset_to_animal = {
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_1.vzg2": "mp16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_3.vzg2": "mp16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_0.vzg2": "mp16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_5.vzg2": "wt_16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_4.vzg2": "wt_16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_2.vzg2": "wt_16_1",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_5.vzg2": "wt_16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_3.vzg2": "wt_16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_4.vzg2": "wt_16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_1.vzg2": "mp16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_2.vzg2": "mp16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_0.vzg2": "mp16_2",
  
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_4.vzg2": "mp11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_3.vzg2": "mp11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_0.vzg2": "mp11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_5.vzg2": "wt11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_1.vzg2": "wt11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_2.vzg2": "wt11_1",
 
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_1.vzg2": "wt11_2",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_5.vzg2": "wt11_2",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_0.vzg2": "wt11_2",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_4.vzg2": "wt11_6",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_3.vzg2": "wt11_6",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_2.vzg2": "wt11_6",
   
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_3.vzg2": "mp11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_4.vzg2": "mp11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_0.vzg2": "wt11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_1.vzg2": "wt11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_2.vzg2": "wt11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_5.vzg2": "mp11_3",
    
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_0.vzg2": "mp8_2",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_1.vzg2": "mp8_1",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_5.vzg2": "mp8_1",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_2.vzg2": "mp8_1",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_3.vzg2": "mp8_2",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_4.vzg2": "mp8_2",        
        
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_1.vzg2": "mp8_3",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_5.vzg2": "mp8_3",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_0.vzg2": "mp8_3",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_3.vzg2": "mp16_4",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_2.vzg2": "mp16_4",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_4.vzg2": "mp16_4",
    
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_5.vzg2": "mp11_4",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_3.vzg2": "mp11_4",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_4.vzg2": "mp11_4",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_2.vzg2": "mp16_3",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_1.vzg2": "mp16_3",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_0.vzg2": "mp16_3"    
}
    
    
# 2) Robust mapping helpers
def _canon(s: str) -> str:
    s = str(s).strip().strip('"').strip("'").replace("\\", "/")
    return s.lower()

def _name(s: str) -> str:
    return PureWindowsPath(str(s)).name.lower()

def _stem(s: str) -> str:
    return PureWindowsPath(str(s)).stem.lower()

def _map_animals_from_datasets(ds_series: pd.Series, mapping: dict) -> pd.Series:
    """Map dataset paths → animal ids via full path, filename, then stem."""
    mapping_full = {_canon(k): v for k, v in mapping.items()}
    mapping_name = {_name(k): v for k, v in mapping.items()}
    mapping_stem = {_stem(k): v for k, v in mapping.items()}

    ds_full = ds_series.astype(str).map(_canon)
    animals = ds_full.map(mapping_full)

    missing = animals.isna()
    if missing.any():
        animals.loc[missing] = ds_full[missing].map(_name).map(mapping_name)

    missing = animals.isna()
    if missing.any():
        animals.loc[missing] = ds_full[missing].map(_stem).map(mapping_stem)

    return animals.fillna("UNKNOWN")

# 3) Build the temporary obs column, plot, save, then clean up
def plot_umap_by_animal_mapping(adata, mapping: dict, datasets_col: str = "Datasets",
                                out_dir: Optional[Path] = None, filename: str = "umap_by_animal_id.pdf",
                                point_size: float = 3.0):
    if datasets_col not in adata.obs:
        # fallback if user earlier mirrored Datasets into section_id
        if "section_id" in adata.obs:
            datasets_col = "section_id"
        else:
            raise KeyError(f"'{datasets_col}' (or 'section_id') not found in adata.obs")

    # Compute temporary labels
    animals = _map_animals_from_datasets(adata.obs[datasets_col], mapping)
    tmp_col = "_animal_id_plot"
    while tmp_col in adata.obs.columns:
        tmp_col += "_x"
    adata.obs[tmp_col] = pd.Categorical(animals)

    # Make sure UMAP exists
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    # Choose output folder
    if out_dir is None:
        # prefer the 'out' Path you already defined; else fallback to OUT_DIR; else CWD
        out_dir = (out if "out" in globals() else OUT_DIR) if "OUT_DIR" in globals() else Path.cwd()
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # Plot and save (figure only)
    ax = sc.pl.umap(
        adata,
        color=[tmp_col],
        legend_loc="right margin",
        size=point_size,
        show=False
    )
    plt.gcf().savefig(out_dir / filename, bbox_inches="tight", dpi=1200)
    plt.close('all')

    # Clean up in-memory only (no files touched)
    del adata.obs[tmp_col]
    # Also drop the auto-added color palette entry if present (purely in-memory)
    adata.uns.pop(f"{tmp_col}_colors", None)

    # Optional console hints
    vc = animals.value_counts().sort_index()
    print("[animal mapping] categories:", list(vc.index))
    print(vc.to_string())

# ---- call it (after UMAP is computed) ----
plot_umap_by_animal_mapping(
    adata,
    dataset_to_animal,
    datasets_col="Datasets",                     # will fall back to 'section_id' if needed
    out_dir=out if "out" in globals() else OUT_DIR,
    filename="umap_by_animal_id.pdf",
    point_size=3.0
)
# =============================================================================


#%% Identify candidate dopaminergic clusters using marker expression

import numpy as np, pandas as pd
from scipy import sparse as sp
import scanpy as sc
import pathlib

adata =  sc.read_h5ad(r"F:/ok/14-10-25_10-07_All_animals_exp_Cells_clustered.h5ad")
# --- genes of interest (case-insensitive lookup) ---
goi = ['Th','Slc6a3', "Ddc", "Nr4a2"]
found = {}
for g in goi:
    hits = [v for v in adata.var_names if str(v).lower()==g.lower()]
    if not hits:
        print(f"WARNING: gene '{g}' not found in adata.var_names")
    else:
        found[g] = hits[0]
genes = list(found.values())
if len(genes) < 2:
    raise ValueError("Need both Th and Slc6a3 present to proceed.")

# --- mean (log1p-normalized) expression per Leiden cluster, only for these genes ---
clusts = list(adata.obs['leiden'].astype('category').cat.categories)
rows = []
for cl in clusts:
    sub = adata[adata.obs['leiden']==cl, genes]
    X = sub.X
    m = np.asarray(X.mean(axis=0)).ravel() if sp.issparse(X) else X.mean(axis=0)
    rows.append(m)

df = pd.DataFrame(rows, index=clusts, columns=genes)
df['sum_two'] = df[genes].sum(1)
df_sorted = df.sort_values('sum_two', ascending=False)
print(df_sorted.head(10))        # top candidate clusters

#saving the dendrogram of expression 
import pathlib, scanpy as sc
out = pathlib.Path(r"F:/vizgen_clustering_out")
out.mkdir(parents=True, exist_ok=True)

# Tell Scanpy where to write figures when using `save=...`
sc.settings.figdir = str(out)

sc.pl.dotplot(
    adata,
    var_names=genes,
    groupby='leiden',
    standard_scale='var',
    dendrogram=True,
    show=False,
    save="_th_slc6a3.pdf"   # writes F:/vizgen_clustering_out/dotplot_th_slc6a3.pdf
)

#%% Export whole-dataset UMAP with fixed point size and highlight the candidate cluster 
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt

adata = sc.read_h5ad(r"F:/ok/14-10-25_10-07_All_animals_exp_Cells_clustered.h5ad")

# choose one point size for ALL UMAPs
pt = 0.2   # or use your _auto_point_size(adata.n_obs)

sc.pl.umap(adata, color="leiden", show=False)
default_palette = adata.uns["leiden_colors"].copy()

# Modify one cluster color
clusters = adata.obs["leiden"].cat.categories
idx = list(clusters).index("14")  # cluster to recolor
default_palette[idx] = "#fe6a97ff" # new red color

fig1 = sc.pl.umap(
    adata,
    color="leiden",
    palette=default_palette,
    s=pt,
    show=False,
    return_fig=True
)

fig1.savefig("F:/vizgen_clustering_out/umap_leiden_dotsizefixed.pdf", bbox_inches="tight")
plt.close(fig1)



# ---------- recoloured "marker_candidates" UMAP ----------
selected = ["14"]
adata.obs["marker_candidates"] = np.where(
    adata.obs["leiden"].isin(selected), "candidate", "other"
)

fig2 = sc.pl.umap(
    adata,
    color="marker_candidates",
    s=pt,   # <-- EXACT SAME VALUE
    palette=["#fe6a97ff", "#D0D0D0"],
    show=False,
    return_fig=True,
)
fig2.savefig("F:/vizgen_clustering_out/candidate_umap_dotsizefixed.pdf", bbox_inches="tight")
plt.close(fig2)

#%% Identify and highlight the candidate dopaminergic parent cluster
import numpy as np
import scanpy as sc

adata = sc.read_h5ad(r"F:/ok/14-10-25_10-07_All_animals_exp_Cells_clustered.h5ad")

# choose your parent cluster(s)
selected = ["14"]   # <- put your Leiden cluster(s) here

# mark candidates vs others (in-memory only)
adata.obs["marker_candidates"] = np.where(
    adata.obs["leiden"].isin(selected), "candidate", "other"
)

# make the UMAP; no saving, just show
fig=sc.pl.umap(
    adata,
    color="marker_candidates",
    s=4,          # or use your _auto_point_size(adata.n_obs)
    palette=["#fe6a97ff", "#9c9c9cff"],  # optional, if you want custom colours
    show=True
)

fig.savefig(
    "candidate_umap.pdf",
    bbox_inches="tight"
)
plt.close(fig)

print("Saved candidate UMAP as candidate_umap.pdf")

#2nd way :
import numpy as np, pandas as pd
from scipy import sparse as sp
import scanpy as sc

adata = sc.read_h5ad(r"F:/ok/14-10-25_10-07_All_animals_exp_Cells_clustered.h5ad")

# helper: case-insensitive gene find
def case_find_genes(adata, names):
    lookup = {str(v).lower(): str(v) for v in adata.var_names}
    found = []
    for g in names:
        k = lookup.get(str(g).lower())
        if k is not None:
            found.append(k)
    return found

# helper: compute marker scores per Leiden cluster
def select_parent_clusters(adata, genes, parent_key="leiden", top_k=1):
    clusts = adata.obs[parent_key].astype("category").cat.categories
    rows = []
    for cl in clusts:
        sub = adata[adata.obs[parent_key] == cl, genes]
        X = sub.X
        m = np.asarray(X.mean(axis=0)).ravel() if sp.issparse(X) else X.mean(axis=0)
        rows.append(m)
    df = pd.DataFrame(rows, index=clusts, columns=genes)
    df["sum_markers"] = df[genes].sum(1)
    selected = df.sort_values("sum_markers", ascending=False).head(top_k).index.tolist()
    return selected

genes = case_find_genes(adata, ["Th", "Slc6a3"])
selected = select_parent_clusters(adata, genes, parent_key="leiden", top_k=1)
print("Selected parent clusters:", selected)

adata.obs["marker_candidates"] = np.where(
    adata.obs["leiden"].isin(selected), "candidate", "other"
)

sc.pl.umap(
    adata,
    color="marker_candidates",
    s=4,
    show=True    # again: just display, no saving
)



#%% Select the dopaminergic parent cluster and perform subclustering

"""
Subclustering pipeline (no gene-expression plots)
- Loads a clustered full AnnData (.h5ad)
- Picks parent clusters using marker means (for selection only)
- Subsets, re-embeds, runs Leiden for the subset
- Saves UMAPs (Windows-safe saving, auto point sizes)
- Writes subset .h5ad
- Optionally writes labels back to the full object and (optionally) saves a Visualizer-ready HDF5

NOTE: This version intentionally omits any gene-expression visualizations (no UMAP-by-gene,
no rank_genes_groups, no marker dotplots, no per-gene CSVs).
"""

# ============================== CONFIG ==============================
# --- I/O ---
INPUT_H5AD = r"F:/ok/14-10-25_10-07_All_animals_exp_Cells_clustered.h5ad"
OUT_DIR    = r"F:/vizgen_clustering_out/subcluster_TH_DAT"

# --- Parent selection (uses expression for selection only; no gene plots) ---
MARKERS    = ("Th", "Slc6a3")  # case-insensitive
TOP_K      = 1                 # choose top-K parent clusters by summed marker means
PARENT_KEY = "leiden"          # parent clustering column in .obs

# --- Embedding/Clustering for the subset ---
N_NEIGHBORS   = 12
N_PCS         = 30
METRIC        = "cosine"
UMAP_MIN_DIST = 0.20
UMAP_SPREAD   = 1.2
UMAP_INIT_POS = "spectral"

# --- Leiden for subset ---
RESOLUTION      = 0.80  # used if TARGET_CLUSTERS is None
TARGET_CLUSTERS = None      # or set e.g. 6 to auto-pick a resolution close to this number

# --- Outputs/behavior ---
WRITE_BACK_TO_FULL       = True   # write subcluster labels back into full AnnData (in-memory)
SAVE_SUBSET_H5AD         = True   # save subset AnnData (.h5ad)
SAVE_FULL_WITH_SUBLABELS = False  # also save the full AnnData with labels
WRITE_VISUALIZER_HDF5    = True   # make a Visualizer-importable HDF5 with the new labels
NEW_LABEL_COL            = "DA_subclusters_Th_Slc6a3"  # column name in .obs for Visualizer

# --- Plotting ---
SHOW_PLOTS = True
SAVE_PLOTS = True
PLOT_DPI   = 500

# --- Misc ---
RANDOM_STATE    = 42
SCALE_MAX_VALUE = 10

# crisp PDFs; points can be rasterized later to keep size small
sc.settings.set_figure_params(dpi=150, dpi_save=1200, vector_friendly=True)

# optional: rasterize the point cloud inside the PDF (huge speed/size win)
RASTERIZE_POINTS = True

# ===================================================================

from pathlib import Path
import re, random
import numpy as np
import pandas as pd
from scipy import sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from pandas.api.types import is_categorical_dtype, is_string_dtype

random.seed(RANDOM_STATE)
np.random.seed(RANDOM_STATE)
sc.settings.autoshow = SHOW_PLOTS
sc.settings.verbosity = 3

# -------------------- plotting helpers (Windows-safe) --------------------
def _safe_fn(name: str) -> str:
    """Keep [A-Za-z0-9_.-], replace others with '_'."""
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name)).strip('_')

def _auto_point_size(n_obs: int) -> float:
    """Reasonable UMAP dot sizes across scales."""
    if n_obs <= 2_000: return 10
    if n_obs <= 10_000: return 6
    if n_obs <= 50_000: return 4
    return max(2, min(3.5, 200000 / max(1, n_obs)))

def _save_current(figpath: Path, dpi=PLOT_DPI):
    figpath.parent.mkdir(parents=True, exist_ok=True)
    plt.gcf().savefig(str(figpath), dpi=dpi, bbox_inches="tight")
    plt.close('all')

# --------------------------- utilities ---------------------------
def _case_find_genes(adata, names):
    """Return list of gene names present in adata.var_names (case-insensitive)."""
    lookup = {str(v).lower(): str(v) for v in adata.var_names}
    found = []
    for g in names:
        k = lookup.get(str(g).lower())
        if k is not None:
            found.append(k)
        else:
            print(f"[WARN] Gene '{g}' not found; skipping.")
    return found

def _select_parent_clusters(adata, genes, parent_key, top_k):
    """Compute mean expression per gene per parent cluster, select top_k by sum."""
    if parent_key not in adata.obs:
        raise KeyError(f"Parent key '{parent_key}' not found in adata.obs")
    clusts = adata.obs[parent_key].astype("category").cat.categories
    rows = []
    for cl in clusts:
        sub = adata[adata.obs[parent_key] == cl, genes]
        X = sub.X
        m = np.asarray(X.mean(axis=0)).ravel() if sp.issparse(X) else X.mean(axis=0)
        rows.append(m)
    df = pd.DataFrame(rows, index=clusts, columns=genes)
    df["sum_markers"] = df[genes].sum(1)
    selected = df.sort_values("sum_markers", ascending=False).head(top_k).index.tolist()
    return df, selected

def _first_present(obs, keys):
    for k in keys:
        if k in obs.columns:
            return k
    return None

def _neighbor_mixing_score(adata, batch):
    """1 - mean(same-batch fraction among neighbors). Higher → better mixing."""
    key = adata.uns.get("neighbors", {}).get("connectivities_key", "connectivities")
    if key not in adata.obsp or batch not in adata.obs:
        return float("nan")
    A = adata.obsp[key].tocsr()
    b = adata.obs[batch].astype("category").cat.codes.to_numpy()
    same = []
    for i in range(A.shape[0]):
        nbrs = A.indices[A.indptr[i]:A.indptr[i+1]]
        if len(nbrs) == 0:
            continue
        same.append((b[nbrs] == b[i]).mean())
    return 1.0 - float(np.mean(same)) if same else float("nan")

def _pick_resolution_for_target(adata, target, grid=(0.80,0.90,1.00,1.05,1.10,1.15,1.20,1.30)):
    """Try a small grid of resolutions; return the one with #clusters closest to target."""
    best = None
    for r in grid:
        sc.tl.leiden(adata, resolution=r, key_added="_tmp_leiden", random_state=RANDOM_STATE)
        k = adata.obs["_tmp_leiden"].nunique()
        score = abs(k - target)
        if best is None or score < best[0]:
            best = (score, r, k)
    adata.obs.drop(columns=["_tmp_leiden"], inplace=True, errors="ignore")
    return best[1], best[2]

def _sanitize_obs_for_h5(df: pd.DataFrame) -> pd.DataFrame:
    """Convert Pandas extension dtypes to plain types for robust HDF5 writing."""
    for col in df.columns:
        if is_string_dtype(df[col]) and not is_categorical_dtype(df[col]):
            df[col] = df[col].astype(str)
        if is_categorical_dtype(df[col]):
            cats = df[col].cat.categories
            if is_string_dtype(cats.dtype):
                df[col] = df[col].astype(object).astype('category')
    return df

# --------------------------- core ---------------------------
def run_subclustering(
    input_path=INPUT_H5AD,
    out_dir=OUT_DIR,
    markers=MARKERS,
    top_k=TOP_K,
    parent_key=PARENT_KEY,
    n_neighbors=N_NEIGHBORS,
    n_pcs=N_PCS,
    metric=METRIC,
    min_dist=UMAP_MIN_DIST,
    spread=UMAP_SPREAD,
    init_pos=UMAP_INIT_POS,
    random_state=RANDOM_STATE,
    resolution=RESOLUTION,
    target_clusters=TARGET_CLUSTERS,
    write_back=WRITE_BACK_TO_FULL,
    save_subset_h5ad=SAVE_SUBSET_H5AD,
    save_full_with_labels=SAVE_FULL_WITH_SUBLABELS,
    write_visualizer_hdf5=WRITE_VISUALIZER_HDF5,
    new_label_col=NEW_LABEL_COL,
):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if SAVE_PLOTS:
        sc.settings.figdir = str(out_dir)

    print(f"[load] {input_path}")
    adata = sc.read_h5ad(str(input_path))

    # convenience: derive section_id from Datasets if present
    if "section_id" not in adata.obs and "Datasets" in adata.obs:
        adata.obs["section_id"] = adata.obs["Datasets"].astype(str)

    # 1) select parent clusters by marker means
    g = _case_find_genes(adata, markers)
    if len(g) == 0:
        raise ValueError("None of the requested marker genes were found in var_names")
    df_scores, selected = _select_parent_clusters(adata, g, parent_key, top_k)
    df_scores.to_csv(out_dir / "parent_cluster_marker_scores.csv")
    print(f"[select] Top {top_k} parent clusters by sum({', '.join(g)}): {selected}")

    # 2) candidate overlay on full UMAP (ONLY if X_umap already exists)
    try:
        if "X_umap" in adata.obsm:
            ad_tmp = adata.copy()
            ad_tmp.obs["marker_candidates"] = np.where(ad_tmp.obs[parent_key].isin(selected), "candidate", "other")
            pt = _auto_point_size(ad_tmp.n_obs)
            sc.pl.umap(ad_tmp, color=["marker_candidates"], s=pt, show=False)
            if SAVE_PLOTS:
                _save_current(out_dir / f"umap_candidates_{_safe_fn('_'.join(g))}.png")
            if SHOW_PLOTS:
                plt.show()
        else:
            print("[info] Full object has no X_umap; skipping candidate overlay.")
    except Exception as e:
        print(f"[warn] candidate UMAP plot skipped: {e}")

    # 3) subset to selected parent clusters
    ad_sub = adata[adata.obs[parent_key].isin(selected)].copy()

    # 4) re-embed subset
    if "X_pca" not in ad_sub.obsm:
        sc.pp.scale(ad_sub, max_value=SCALE_MAX_VALUE)
        sc.tl.pca(ad_sub, n_comps=n_pcs, random_state=random_state)
    sc.pp.neighbors(ad_sub, n_neighbors=n_neighbors, n_pcs=n_pcs, metric=metric)
    sc.tl.umap(ad_sub, min_dist=min_dist, spread=spread, init_pos=init_pos, random_state=random_state)

    # 5) subset Leiden
    if target_clusters is not None:
        best_res, achieved = _pick_resolution_for_target(ad_sub, target_clusters)
        print(f"[leiden] auto-picked resolution={best_res} → {achieved} clusters (target {target_clusters})")
        sc.tl.leiden(ad_sub, resolution=best_res, key_added="leiden_sub", random_state=random_state)
    else:
        sc.tl.leiden(ad_sub, resolution=resolution, key_added="leiden_sub", random_state=random_state)

   # 6) UMAP of subclusters (only leiden_sub) — save as true PDF
    pt = _auto_point_size(ad_sub.n_obs)
    
    # build the figure but do NOT show it yet
    fig = sc.pl.umap(
        ad_sub,
        color="leiden_sub",
        s=pt,
        show=False,
        return_fig=True,   # <-- get the figure object to save reliably
    )
    
    # (optional) rasterize the dot cloud inside the PDF to reduce size
    if 'RASTERIZE_POINTS' in globals() and RASTERIZE_POINTS:
        for ax in fig.axes:
            for coll in getattr(ax, "collections", []):
                try:
                    coll.set_rasterized(True)
                except Exception:
                    pass
    
    # save as PDF explicitly (don’t rely on Scanpy’s save=…)
    pdf_path = Path(out_dir) / "umap_leiden_sub.pdf"
    fig.savefig(pdf_path, bbox_inches="tight")   # dpi from set_figure_params is used for raster parts
    print(f"[save] {pdf_path}")
    
    # optional on-screen preview
    if SHOW_PLOTS:
        plt.show()
    plt.close(fig)

    # 7) diagnostics (sizes, silhouette, batch mixing)
    sizes = pd.Series(ad_sub.obs["leiden_sub"]).value_counts().sort_index()
    sizes.to_csv(out_dir / "cluster_sizes.csv")

    if "X_pca" in ad_sub.obsm:
        from sklearn.metrics import silhouette_score
        labs = ad_sub.obs["leiden_sub"].astype("category").cat.codes.to_numpy()
        sil = silhouette_score(ad_sub.obsm["X_pca"], labs)
        with open(out_dir / "silhouette.txt", "w") as f:
            f.write(f"silhouette_on_pca: {sil:.4f}\n")
        print(f"[silhouette] {sil:.4f}")

    bk = _first_present(ad_sub.obs, ("section_id", "Datasets", "slide_id", "animal_id"))
    if bk:
        mix = _neighbor_mixing_score(ad_sub, bk)
        with open(out_dir / "batch_diagnostics.txt", "w") as f:
            f.write(f"neighbor_mixing({bk}) = {mix}\n")
        print(f"[batch] neighbor_mixing({bk}) = {mix:.3f}")

    # Ensure log1p metadata exists (some Scanpy versions expect this)
    if 'log1p' not in ad_sub.uns or 'base' not in ad_sub.uns.get('log1p', {}):
        ad_sub.uns['log1p'] = {'base': None}

    # 8) write subset object
    if save_subset_h5ad:
        sub_path = Path(out_dir) / ("subset_" + "_".join(g) + ".h5ad")
        ad_sub.write_h5ad(sub_path)
        print(f"[write] subset h5ad → {sub_path}")

    # 9) optionally write subcluster labels back to full object & (optionally) save
    if write_back:
        col = "leiden_sub_" + "_".join(g)
        adata.obs[col] = pd.NA
        adata.obs.loc[ad_sub.obs_names, col] = ad_sub.obs["leiden_sub"].astype(str).values
        adata.obs["leiden_parent_sub"] = adata.obs[parent_key].astype(str)
        mask = adata.obs[col].notna()
        adata.obs.loc[mask, "leiden_parent_sub"] = (
            adata.obs.loc[mask, parent_key].astype(str)
                .str.cat(adata.obs.loc[mask, col].astype(str), sep="_")
        )
        if save_full_with_labels:
            out_full = Path(out_dir) / (Path(input_path).stem + "_with_subclusters.h5ad")
            adata.write_h5ad(out_full)
            print(f"[write] full h5ad with labels → {out_full}")

    # 10) optional: Visualizer-ready HDF5 of the full object with the new labels
    if write_visualizer_hdf5:
        cats = None
        if new_label_col not in adata.obs:
            # mirror from our main sublabel if user kept write_back=True
            # otherwise create a copy using the same logic
            col = "leiden_sub_" + "_".join(g)
            if col in adata.obs:
                adata.obs[new_label_col] = adata.obs[col].astype("category")
            else:
                labels_full = pd.Series(index=adata.obs_names, dtype=object)
                labels_sub = ad_sub.obs["leiden_sub"].astype(str)
                labels_full.loc[labels_sub.index.intersection(adata.obs_names)] = labels_sub.values
                labels_full = labels_full.fillna("Unclustered")
                adata.obs[new_label_col] = pd.Categorical(labels_full)
        cats = list(adata.obs[new_label_col].cat.categories)

        # colors
        cmap = plt.get_cmap("tab20")
        adata.uns[f"{new_label_col}_colors"] = [to_hex(cmap(i % cmap.N)) for i in range(len(cats))]

        # sanitize obs for HDF5
        adata.obs = _sanitize_obs_for_h5(adata.obs)

        out_h5 = Path(out_dir) / "visualizer_import_subclusters.hdf5"
        # uncompressed write (often safer for Visualizer)
        adata.write_h5ad(out_h5)
        print(f"[write] Visualizer HDF5 → {out_h5}")

    print(
        f"[settings] n_neighbors={n_neighbors}, n_pcs={n_pcs}, "
        f"min_dist={min_dist}, spread={spread}, metric='{metric}', init='{init_pos}', seed={random_state}"
    )
    print("[clusters] n_leiden_sub =", ad_sub.obs["leiden_sub"].nunique())
    print("[done] Subclustering complete.")

    return ad_sub, selected

# --------------------------- run ---------------------------
if __name__ == "__main__":
    adata_sub, selected_clusters = run_subclustering()
    print("Selected parent clusters:", selected_clusters)


#%% Plot the five genes with highest mean expression per subcluster

# ==== Top-5 genes PER subcluster — barplot grid (exactly 5 per cluster) ====
from pathlib import Path
from typing import Optional, Tuple, Sequence
import numpy as np, pandas as pd, scanpy as sc, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy import sparse as sp
import math

def gray_bordeaux_cmap():
    return LinearSegmentedColormap.from_list(
        "gray_bordeaux", ["#E6E6E6", "#BDBDBD", "#800020", "#2B000F"], N=256
    )

def _pick_groupby(adata, prefer: Optional[str] = None) -> str:
    if prefer and prefer in adata.obs: return prefer
    if "leiden_sub" in adata.obs:     return "leiden_sub"
    if "leiden" in adata.obs:         return "leiden"
    raise KeyError("No 'leiden_sub' or 'leiden' found in .obs — pass groupby=...")

def _choose_source(adata, layer: Optional[str]):
    # returns (use_layer, use_raw)
    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not in adata.layers")
        return layer, False
    if getattr(adata, "raw", None) is not None and adata.raw is not None:
        return None, True
    return None, False

def _get_matrix_and_varnames(adata, mask, use_layer, use_raw):
    if use_raw:
        return adata.raw[mask, :].X, adata.raw.var_names
    if use_layer is not None:
        return adata.layers[use_layer][mask, :], adata.var_names
    return adata.X[mask, :], adata.var_names

def _mean_expr(X):
    if sp.issparse(X):
        return np.asarray(X.mean(axis=0)).ravel()
    return np.asarray(X).mean(axis=0)

# --- replace your _exclude_mask with this (faster + robust) ---
def _exclude_mask(varnames, keywords):
    """Boolean mask for names containing any keyword (case-insensitive)."""
    if not keywords:
        return np.zeros(len(varnames), dtype=bool)
    names = np.asarray(varnames, dtype=str)        # ensure plain strings
    names_lower = np.char.lower(names)
    mask = np.zeros(names_lower.size, dtype=bool)
    for kw in keywords:
        if kw:
            k = str(kw).lower()
            # simple & robust: list comprehension, no np.char ops on object arrays
            mask |= np.array([k in s for s in names_lower], dtype=bool)
    return mask

def plot_top5_bars_per_cluster(
    adata,
    groupby: Optional[str] = None,
    layer: Optional[str] = None,
    exclude_keywords: Tuple[str, ...] = ("Blank",),
    exclude_genes: Tuple[str, ...] = ("Th", "Slc6a3"),  # exact-name excludes
    out_pdf: str = r"F:/vizgen_clustering_out/subcluster_TH_DAT/top5_per_subcluster_bars.pdf",
    share_x: bool = False,
    figsize_per_col: float = 3.6,
    figsize_per_row: float = 2.8,
    max_cols: int = 4,
    title_prefix: str = "Top-5 (mean expr)"
) -> Path:
    gb = _pick_groupby(adata, prefer=groupby)
    use_layer, use_raw = _choose_source(adata, layer)
    print(f"[info] groupby='{gb}' | source={'raw' if use_raw else ('layer:'+use_layer if use_layer else 'X')}")

    cats = adata.obs[gb].astype("category").cat.categories
    per_cluster, global_max = [], 0.0

    ex_names_set = {str(x).lower() for x in (exclude_genes or ())}

    for c in cats:
        mask = (adata.obs[gb] == c).to_numpy()
        if not mask.any():
            per_cluster.append((c, [], []))
            continue

        Xc, vnames_idx = _get_matrix_and_varnames(adata, mask, use_layer, use_raw)
        vnames = np.asarray(vnames_idx, dtype=str)       # <-- ensure plain string dtype

        means = _mean_expr(Xc)

        # 1) keyword-based exclusion
        ex_kw = _exclude_mask(vnames, exclude_keywords)

        # 2) exact-name exclusion
        v_lower = np.char.lower(vnames)
        ex_exact = np.isin(v_lower, list(ex_names_set)) if ex_names_set else np.zeros_like(ex_kw)

        ex_all = ex_kw | ex_exact
        means = means.copy()
        means[ex_all] = -np.inf

        # pick EXACTLY 5 where possible
        valid_idx = np.where(np.isfinite(means))[0]
        if valid_idx.size == 0:
            per_cluster.append((c, [], []))
            continue

        top_k = 5
        if valid_idx.size <= top_k:
            top_idx = valid_idx[np.argsort(means[valid_idx])[::-1]]
            if valid_idx.size < top_k:
                print(f"[warn] cluster {c}: only {valid_idx.size} usable genes after exclusions (showing all).")
        else:
            part = np.argpartition(means[valid_idx], kth=valid_idx.size - top_k)[-top_k:]
            top_idx = valid_idx[part[np.argsort(means[valid_idx][part])[::-1]]]

        genes = [vnames[i] for i in top_idx]
        vals  = [float(means[i]) for i in top_idx]

        per_cluster.append((c, genes, vals))
        if vals:
            global_max = max(global_max, max(vals))

    # ---- figure grid (unchanged) ----
    n_clusters = len(cats)
    ncols = min(max_cols, n_clusters if n_clusters > 0 else 1)
    nrows = int(math.ceil(n_clusters / ncols))
    fig_w = max(6.0, ncols * figsize_per_col)
    fig_h = max(4.0, nrows * figsize_per_row)
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h), squeeze=False)
    cmap = gray_bordeaux_cmap()

    for idx_c, (c, genes, vals) in enumerate(per_cluster):
        r, k = divmod(idx_c, ncols)
        ax = axes[r, k]
        ax.set_title(f"{title_prefix}\ncluster {c} (n={int((adata.obs[gb]==c).sum())})", fontsize=10)

        if not genes:
            ax.text(0.5, 0.5, "no data", ha="center", va="center", fontsize=10)
            ax.set_axis_off()
            continue

        order = np.argsort(vals)[::-1]
        genes_ord = [genes[i] for i in order]
        vals_ord  = [vals[i] for i in order]

        vmax = (global_max if share_x else max(vals_ord))
        vmin = 0.0
        norm = Normalize(vmin=vmin, vmax=vmax if vmax > vmin else vmin + 1e-9)
        colors = [cmap(norm(v)) for v in vals_ord]

        ax.barh(genes_ord, vals_ord, color=colors, edgecolor="black", linewidth=0.5)
        ax.invert_yaxis()
        ax.set_xlabel("mean expression")
        if share_x:
            ax.set_xlim(vmin, vmax * 1.05 if vmax > 0 else 1.0)
        ax.tick_params(axis='y', labelsize=9)
        ax.tick_params(axis='x', labelsize=8)

    for idx_c in range(n_clusters, nrows * ncols):
        r, k = divmod(idx_c, ncols)
        axes[r, k].set_axis_off()

    fig.tight_layout()
    out_pdf = Path(out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_pdf), bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_pdf}")
    return out_pdf

pdf = plot_top5_bars_per_cluster(
    adata_sub,
    groupby=None,
    layer=None,
    exclude_keywords=("Blank",),
    exclude_genes=("Th","Slc6a3"),
    out_pdf=r"F:/vizgen_clustering_out/subcluster_TH_DAT/top5_per_subcluster_bars.pdf",
    share_x=False,
    max_cols=4
)
print("Saved:", pdf)

#%% Optional: Recover expression layers and compare subcluster gene rankings

import scanpy as sc
from scipy import sparse as sp
import numpy as np
import pandas as pd

RAW_H5 = r"F:/25-07-25_10-04_all_mitopark_Cells.hdf5"  # your original export

def attach_layers_for_plots(adata_sub, raw_h5):
    ad_raw = sc.read_h5ad(raw_h5)

    # 1) align cells (subset must be present in raw)
    missing_cells = set(adata_sub.obs_names) - set(ad_raw.obs_names)
    if missing_cells:
        raise ValueError(f"{len(missing_cells)} cells from subset are missing in RAW file.")

    ad_raw = ad_raw[adata_sub.obs_names, :].copy()

    # 2) align genes (we expect the same panel; warn if not)
    miss_genes = set(adata_sub.var_names) - set(ad_raw.var_names)
    if miss_genes:
        # If this ever happens, we’ll just drop those genes from plotting
        # but in most MERSCOPE workflows this should be empty.
        print(f"[warn] RAW file missing {len(miss_genes)} genes used in subset (they won't have counts/log1p).")

    # Reorder RAW to the subset’s gene order, keep only common genes
    common = adata_sub.var_names.intersection(ad_raw.var_names)
    # maps: subset gene order indices for 'common'
    take_cols = pd.Index(ad_raw.var_names).get_indexer(common)
    # build a full-size counts matrix with zeros, then fill common columns
    n, p = adata_sub.n_obs, adata_sub.n_vars
    # choose dtype from raw
    dtype = ad_raw.X.dtype if hasattr(ad_raw.X, "dtype") else np.float32

    # Start from zeros; fill only columns that actually exist in RAW
    counts_full = sp.csr_matrix((n, p), dtype=dtype)
    if len(common) > 0:
        # Extract counts for common genes in subset cell order, as dense or sparse
        raw_common = ad_raw[:, common].X
        # place them into the correct columns of counts_full
        # convert to CSR for safe hstack
        if not sp.issparse(raw_common):
            raw_common = sp.csr_matrix(raw_common)
        # column positions in the subset corresponding to 'common'
        target_pos = pd.Index(adata_sub.var_names).get_indexer(common)
        # Build as hstack of blocks is cumbersome; a simple approach:
        # create a LIL (efficient assignment), then convert back
        counts_full = sp.lil_matrix((n, p), dtype=raw_common.dtype)
        counts_full[:, target_pos] = raw_common
        counts_full = counts_full.tocsr()

    # Save layers (counts)
    adata_sub.layers["counts"] = counts_full

    # 3) log1p on normalized counts (unscaled)
    # work on a tiny AnnData wrapper to reuse scanpy normalize/log1p safely
    ad_tmp = sc.AnnData(counts_full.copy(), var=adata_sub.var.copy(), obs=adata_sub.obs.copy())
    sc.pp.normalize_total(ad_tmp)
    sc.pp.log1p(ad_tmp)
    adata_sub.layers["log1p"] = ad_tmp.X.copy()

    print("Layers attached (in-memory only):", list(adata_sub.layers.keys()))
    return adata_sub

# Run it once
attach_layers_for_plots(adata_sub, RAW_H5)

# 1) Best for visuals: unscaled log1p
pdf1 = plot_top5_bars_per_cluster(
    adata_sub,
    layer="log1p",
    exclude_keywords=("Blank",),
    exclude_genes=("Th","Slc6a3"),
    out_pdf=r"F:/vizgen_clustering_out/subcluster_TH_DAT/top5__log1p.pdf",
    title_prefix="Top-5 (mean log1p, norm)"
)

# 2) Raw counts (sanity/model checks)
pdf2 = plot_top5_bars_per_cluster(
    adata_sub,
    layer="counts",
    exclude_keywords=("Blank",),
    exclude_genes=("Th","Slc6a3"),
    out_pdf=r"F:/vizgen_clustering_out/subcluster_TH_DAT/top5__counts.pdf",
    title_prefix="Top-5 (mean counts)"
)

# 3) Your current scaled matrix (.X) = z-scored log1p
pdf3 = plot_top5_bars_per_cluster(
    adata_sub,
    layer=None,  # falls back to .X
    exclude_keywords=("Blank",),
    exclude_genes=("Th","Slc6a3"),
    out_pdf=r"F:/vizgen_clustering_out/subcluster_TH_DAT/top5__scaledX.pdf",
    title_prefix="Top-5 (mean z-score of log1p)"
)

print("Saved:\n", pdf1, "\n", pdf2, "\n", pdf3)


#%% Export full-dataset and subcluster UMAPs for Vizgen Vizualizer

from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from pandas.api.types import is_categorical_dtype, is_string_dtype, is_integer_dtype, is_bool_dtype

# ---- paths (edit if needed) ----
FULL_H5AD   = r"F:/ok/14-10-25_10-07_All_animals_exp_Cells_clustered.h5ad"
SUBSET_H5AD = r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad"

OUT_FULL_H5   = r"F:/vizgen_clustering_out/subcluster_TH_DAT/visualizer_full_with_DA_subclusters_viz_ready.hdf5"
OUT_SUBSET_H5 = r"F:/vizgen_clustering_out/subcluster_TH_DAT/visualizer_subset_DA_subclusters_viz_ready.hdf5"

# the subcluster column name(s) that may exist in the subset
CANDIDATE_SUB_LABELS = ["leiden_sub_Th_Slc6a3", "leiden_parent_sub", "leiden_sub", "leiden"]
# how it should appear in the FULL file
NEW_LABEL_COL = "DA_subclusters_Th_Slc6a3"

# ---- helpers ----
def _ensure_cell_id(ad):
    for k in ["Cell ID","cell_id","Cell_ID","cellID","CellId","cellId","ID","cell_index"]:
        if k in ad.obs.columns:
            ad.obs[k] = ad.obs[k].astype(str)
            return k
    ad.obs["Cell ID"] = ad.obs_names.astype(str)
    return "Cell ID"

def _sanitize_obs(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        s = out[c]
        # extension dtypes → plain str to avoid importer issues
        if str(s.dtype) in ("Int64","UInt64","UInt32","UInt16","UInt8","boolean"):
            out[c] = s.astype(str)
        elif is_string_dtype(s) and not is_categorical_dtype(s):
            out[c] = s.astype(str)
        elif is_categorical_dtype(s):
            cats = s.cat.categories
            if is_string_dtype(cats.dtype):
                out[c] = s.astype(object).astype("category")
    return out

def _sanitize_var(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        s = out[c]
        if is_string_dtype(s) and not is_categorical_dtype(s):
            out[c] = s.astype(str)
        elif is_integer_dtype(s) and str(s.dtype) != "int64":
            out[c] = s.astype("int64", copy=False)
        elif is_bool_dtype(s):
            out[c] = s.astype(bool, copy=False)
    return out

def _ensure_colors(ad, label_key: str):
    cats = list(ad.obs[label_key].astype("category").cat.categories)
    cmap = plt.get_cmap("tab20")
    ad.uns[f"{label_key}_colors"] = [to_hex(cmap(i % cmap.N)) for i in range(len(cats))]
    return cats

def _strip_heavy_bits(ad):
    # Drop things some Visualizer builds choke on
    ad.uns.pop("neighbors", None)
    for k in ("connectivities", "distances"):
        if k in ad.obsp:
            del ad.obsp[k]
    try:
        ad.layers.clear()
    except Exception:
        for k in list(ad.layers.keys()):
            del ad.layers[k]
    try:
        ad.raw = None
    except Exception:
        pass

# ---- create FULL viz-ready: global UMAP + subcluster labels on all cells ----
def export_full_with_sublabels(full_h5ad: str, subset_h5ad: str, out_h5: str, new_label_col: str):
    ad_full = sc.read_h5ad(full_h5ad)
    ad_full.var_names_make_unique()

    # find sublabel column on the subset
    ad_sub = sc.read_h5ad(subset_h5ad)
    sub_label = next((c for c in CANDIDATE_SUB_LABELS if c in ad_sub.obs.columns), None)
    if sub_label is None:
        raise KeyError(f"No subcluster label column found in subset. Have: {list(ad_sub.obs.columns)[:30]}")

    # ensure Cell ID
    _ensure_cell_id(ad_full)

    # build full-length label Series (non-subset → 'Unclustered')
    s_sub = ad_sub.obs[sub_label].astype(str)
    s_full = pd.Series("Unclustered", index=ad_full.obs_names, dtype=object)
    idx = s_full.index.intersection(s_sub.index)
    s_full.loc[idx] = s_sub.loc[idx].values
    ad_full.obs[new_label_col] = pd.Categorical(s_full)

    # mirror to 'Custom Cells groups' (many Visualizer builds look for this)
    ad_full.obs["Custom Cells groups"] = ad_full.obs[new_label_col].copy()

    # ensure colors for both keys
    _ensure_colors(ad_full, new_label_col)
    _ensure_colors(ad_full, "Custom Cells groups")

    # UMAP must be present; cast to float64 Nx2
    if "X_umap" not in ad_full.obsm:
        raise KeyError("Full object has no obsm['X_umap']; compute UMAP before exporting.")
    ad_full.obsm["X_umap"] = np.asarray(ad_full.obsm["X_umap"])[:, :2].astype(np.float64, copy=False)

    # sanitize & strip
    ad_full.obs = _sanitize_obs(ad_full.obs)
    ad_full.var = _sanitize_var(ad_full.var)
    _strip_heavy_bits(ad_full)

    # write UNCOMPRESSED
    Path(out_h5).parent.mkdir(parents=True, exist_ok=True)
    ad_full.write_h5ad(out_h5)
    print(f"[write] FULL viz-ready → {out_h5}  (label='{new_label_col}', categories={ad_full.obs[new_label_col].cat.categories.size})")

# ---- create SUBSET viz-ready: subcluster UMAP + subcluster labels ----
from typing import Optional
import numpy as np, pandas as pd, scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
from pandas.api.types import is_categorical_dtype, is_string_dtype, is_integer_dtype, is_bool_dtype
from pathlib import Path

def export_subset_subumap(subset_h5ad: str, out_h5: str, prefer_label: Optional[str] = None):
    """
    Write a *subset-only* viz-ready HDF5 whose embedding is the SUBCLUSTER UMAP.
    Ensures:
      - obsm['X_umap'] is float64, shape (n,2)
      - obs['umap_1'], obs['umap_2'] exist (floats)
      - obs['leiden'] holds the subcluster labels (categorical strings)
      - obs['Custom Cells groups'] mirrors 'leiden'
      - uns['leiden_colors'] and uns['Custom Cells groups_colors'] exist with matching lengths
      - 'Cell ID' column exists (plain strings)
      - heavy bits removed; uncompressed write
    """
    # ---- load subset ----
    ad_sub = sc.read_h5ad(subset_h5ad)
    ad_sub.var_names_make_unique()

    # pick the label column
    CANDIDATE_SUB_LABELS = ["leiden_sub_Th_Slc6a3", "leiden_parent_sub", "leiden_sub", "leiden"]
    label_key = prefer_label if (prefer_label and prefer_label in ad_sub.obs) else \
                next((c for c in CANDIDATE_SUB_LABELS if c in ad_sub.obs.columns), None)
    if label_key is None:
        raise KeyError(f"No subcluster label column found in subset. Have: {list(ad_sub.obs.columns)[:30]}")

    # make 'leiden' be *the* subcluster labels (categorical strings)
    ad_sub.obs["leiden"] = ad_sub.obs[label_key].astype(str).astype("category")
    ad_sub.obs["Custom Cells groups"] = ad_sub.obs["leiden"].copy()

    # colors for both keys
    def _ensure_colors(ad, key):
        cats = list(ad.obs[key].cat.categories)
        cmap = plt.get_cmap("tab20")
        ad.uns[f"{key}_colors"] = [to_hex(cmap(i % cmap.N)) for i in range(len(cats))]
    _ensure_colors(ad_sub, "leiden")
    _ensure_colors(ad_sub, "Custom Cells groups")

    # ensure Cell ID (plain strings)
    def _ensure_cell_id(ad):
        for k in ["Cell ID","cell_id","Cell_ID","cellID","CellId","cellId","ID","cell_index"]:
            if k in ad.obs.columns:
                ad.obs[k] = ad.obs[k].astype(str)
                return k
        ad.obs["Cell ID"] = ad.obs_names.astype(str)
        return "Cell ID"
    _ensure_cell_id(ad_sub)

    # embedding: must exist and be float64 (n×2); also copy into obs columns
    if "X_umap" not in ad_sub.obsm:
        raise KeyError("Subset object has no obsm['X_umap']; compute subcluster UMAP before exporting.")
    um = np.asarray(ad_sub.obsm["X_umap"])[:, :2].astype(np.float64, copy=False)
    ad_sub.obsm["X_umap"] = um
    ad_sub.obs["umap_1"] = um[:, 0]
    ad_sub.obs["umap_2"] = um[:, 1]
    for k in list(ad_sub.obsm.keys()):
        if k != "X_umap":
            del ad_sub.obsm[k]

    # sanitize obs/var dtypes for HDF5 readers
    from pandas.api.types import is_categorical_dtype, is_string_dtype, is_integer_dtype, is_bool_dtype
    def _sanitize_obs(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        for c in out.columns:
            s = out[c]
            if str(s.dtype) in ("Int64","UInt64","UInt32","UInt16","UInt8","boolean"):
                out[c] = s.astype(str)
            elif is_string_dtype(s) and not is_categorical_dtype(s):
                out[c] = s.astype(str)
            elif is_categorical_dtype(s):
                cats = s.cat.categories
                if is_string_dtype(cats.dtype):
                    out[c] = s.astype(object).astype("category")
        return out
    def _sanitize_var(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        for c in out.columns:
            s = out[c]
            if is_string_dtype(s) and not is_categorical_dtype(s):
                out[c] = s.astype(str)
            elif is_integer_dtype(s) and str(s.dtype) != "int64":
                out[c] = s.astype("int64", copy=False)
            elif is_bool_dtype(s):
                out[c] = s.astype(bool, copy=False)
        return out
    ad_sub.obs = _sanitize_obs(ad_sub.obs)
    ad_sub.var = _sanitize_var(ad_sub.var)

    # drop heavy bits that sometimes crash Visualizer
    ad_sub.uns.pop("neighbors", None)
    for k in ("connectivities", "distances"):
        if k in ad_sub.obsp:
            del ad_sub.obsp[k]
    try:
        ad_sub.layers.clear()
    except Exception:
        for k in list(ad_sub.layers.keys()):
            del ad_sub.layers[k]
    try:
        ad_sub.raw = None
    except Exception:
        pass

    # write UNCOMPRESSED
    out_h5 = str(out_h5)
    Path(out_h5).parent.mkdir(parents=True, exist_ok=True)
    ad_sub.write_h5ad(out_h5)
    print(f"[write] SUBSET viz-ready → {out_h5}  (cells={ad_sub.n_obs}, k={ad_sub.obs['leiden'].cat.categories.size})")

# --- actually run the exports ---
sc.settings.verbosity = 3  # optional: more logging

# 1) FULL file: global UMAP + subcluster labels painted across all cells
export_full_with_sublabels(
    full_h5ad=FULL_H5AD,
    subset_h5ad=SUBSET_H5AD,
    out_h5=OUT_FULL_H5,
    new_label_col=NEW_LABEL_COL
)

# 2) SUBSET file: only the selected cells, with the SUBCLUSTER UMAP
# If you want to force a specific sublabel column, pass prefer_label="leiden_sub_Th_Slc6a3"
export_subset_subumap(
    subset_h5ad=SUBSET_H5AD,
    out_h5=OUT_SUBSET_H5,
    prefer_label=None  # or "leiden_sub_Th_Slc6a3"
)

# --- quick sanity check (optional) ---
ad_full_chk = sc.read_h5ad(OUT_FULL_H5)
print("[FULL] cells:", ad_full_chk.n_obs, "| X_umap:", ad_full_chk.obsm["X_umap"].shape,
      "| label col:", NEW_LABEL_COL, " (", ad_full_chk.obs[NEW_LABEL_COL].cat.categories.size, "cats )")

ad_sub_chk = sc.read_h5ad(OUT_SUBSET_H5)
print("[SUBSET] cells:", ad_sub_chk.n_obs, "| X_umap:", ad_sub_chk.obsm["X_umap"].shape,
      "| leiden cats:", ad_sub_chk.obs["leiden"].cat.categories.size)
print("subset umap head:\n", np.round(ad_sub_chk.obsm["X_umap"][:3], 3))


#%% Inspect animal distribution across the dopaminergic subcluster UMAP
# ==== UMAP colored by dataset→animal mapping (SUBCLUSTER view) ====
from typing import Union
from pathlib import PureWindowsPath, Path
import scanpy as sc, matplotlib.pyplot as plt
import pandas as pd, numpy as np

# make sure figures save even if autoshow=True elsewhere
_sc_autoshow_prev = sc.settings.autoshow
sc.settings.autoshow = False
sc.settings.set_figure_params(dpi=150, dpi_save=1200, vector_friendly=True)

RASTERIZE_POINTS = True  # smaller/faster PDFs

def _canon(s: str) -> str:
    s = str(s).strip().strip('"').strip("'").replace("\\", "/")
    return s.lower()

def _name(s: str) -> str:
    return PureWindowsPath(str(s)).name.lower()

def _stem(s: str) -> str:
    return PureWindowsPath(str(s)).stem.lower()

def _map_animals_from_datasets(ds_series: pd.Series, mapping: dict) -> pd.Series:
    mapping_full = {_canon(k): v for k, v in mapping.items()}
    mapping_name = {_name(k): v for k, v in mapping.items()}
    mapping_stem = {_stem(k): v for k, v in mapping.items()}
    ds_full = ds_series.astype(str).map(_canon)
    animals = ds_full.map(mapping_full)
    miss = animals.isna()
    if miss.any():
        animals.loc[miss] = ds_full[miss].map(_name).map(mapping_name)
    miss = animals.isna()
    if miss.any():
        animals.loc[miss] = ds_full[miss].map(_stem).map(mapping_stem)
    return animals.fillna("UNKNOWN")

def _auto_point_size(n_obs: int) -> float:
    if n_obs <= 2_000: return 10
    if n_obs <= 10_000: return 6
    if n_obs <= 50_000: return 4
    return max(2, min(3.5, 200000 / max(1, n_obs)))

def plot_sub_umap_by_animal_mapping(
    adata, mapping: dict, datasets_col: str = "Datasets",
    out_dir: Union[Path, str] = None, filename: str = "sub_umap_by_animal_id.pdf",
    point_size: float = None
) -> Path:
    # pick column
    if datasets_col not in adata.obs:
        if "section_id" in adata.obs:
            datasets_col = "section_id"
        else:
            raise KeyError(f"'{datasets_col}' (or 'section_id') not found in adata.obs")

    # temp categorical + colors
    animals = _map_animals_from_datasets(adata.obs[datasets_col], mapping)
    tmp_col = "__tmp_animal_id__"
    while tmp_col in adata.obs.columns:
        tmp_col += "_x"
    cats = sorted(animals.unique().tolist())
    adata.obs[tmp_col] = pd.Categorical(animals, categories=cats)
    cmap = plt.get_cmap("tab20")
    adata.uns[f"{tmp_col}_colors"] = [plt.matplotlib.colors.to_hex(cmap(i % cmap.N)) for i in range(len(cats))]

    # ensure UMAP
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    # output path
    if out_dir is None:
        out_dir = OUT_DIR if "OUT_DIR" in globals() else Path("F:/vizgen_clustering_out/subcluster_TH_DAT")
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / filename

    # plot → figure object → save
    s = _auto_point_size(adata.n_obs) if point_size is None else float(point_size)
    obj = sc.pl.umap(adata, color=tmp_col, legend_loc="right margin", s=s, show=False, return_fig=True)

    # handle both possible return types (Figure or Axes)
    fig = obj if hasattr(obj, "savefig") else plt.gcf()

    if RASTERIZE_POINTS:
        for ax in fig.axes:
            for coll in getattr(ax, "collections", []):
                try: coll.set_rasterized(True)
                except Exception: pass

    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {pdf_path}")

    # cleanup (in-memory only)
    del adata.obs[tmp_col]
    adata.uns.pop(f"{tmp_col}_colors", None)

    # quick summary
    vc = animals.value_counts().sort_index()
    print("[animal mapping] categories:", list(vc.index))
    print(vc.to_string())
    return pdf_path

# restore previous autoshow behavior
sc.settings.autoshow = _sc_autoshow_prev

# 1) your mapping dict must be defined
dataset_to_animal = {
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_1.vzg2": "mp16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_3.vzg2": "mp16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_0.vzg2": "mp16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_5.vzg2": "wt_16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_4.vzg2": "wt_16_1",
    "F:/all_vizgen_files/202507171430_Mitopark16w2ko_VMSC13902_region_2.vzg2": "wt_16_1",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_5.vzg2": "wt_16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_3.vzg2": "wt_16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_4.vzg2": "wt_16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_1.vzg2": "mp16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_2.vzg2": "mp16_2",
    "F:/all_vizgen_files/202507191504_Mitopark16kowt_VMSC13902_region_0.vzg2": "mp16_2",
  
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_4.vzg2": "mp11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_3.vzg2": "mp11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_0.vzg2": "mp11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_5.vzg2": "wt11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_1.vzg2": "wt11_1",
    "F:/all_vizgen_files/202507021142_MITOPARK1_VMSC13902_region_2.vzg2": "wt11_1",
 
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_1.vzg2": "wt11_2",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_5.vzg2": "wt11_2",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_0.vzg2": "wt11_2",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_4.vzg2": "wt11_6",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_3.vzg2": "wt11_6",
    "F:/all_vizgen_files/202507101144_Mitopark2_VMSC13902_region_2.vzg2": "wt11_6",
   
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_3.vzg2": "mp11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_4.vzg2": "mp11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_0.vzg2": "wt11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_1.vzg2": "wt11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_2.vzg2": "wt11_3",
    "F:/all_vizgen_files/202507121804_mitopark3_VMSC13902_region_5.vzg2": "mp11_3",
    
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_0.vzg2": "mp8_2",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_1.vzg2": "mp8_1",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_5.vzg2": "mp8_1",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_2.vzg2": "mp8_1",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_3.vzg2": "mp8_2",
    "F:/all_vizgen_files/202509291034_DMC-Mitopark-8weeks_VMSC13902_region_4.vzg2": "mp8_2",        
        
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_1.vzg2": "mp8_3",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_5.vzg2": "mp8_3",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_0.vzg2": "mp8_3",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_3.vzg2": "mp16_4",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_2.vzg2": "mp16_4",
    "F:/all_vizgen_files/202510021003_dmc-mitopark-8wks16wks_VMSC13902_region_4.vzg2": "mp16_4",
    
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_5.vzg2": "mp11_4",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_3.vzg2": "mp11_4",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_4.vzg2": "mp11_4",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_2.vzg2": "mp16_3",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_1.vzg2": "mp16_3",
    "F:/all_vizgen_files/202509261129_mitopark10-16_VMSC13902_region_0.vzg2": "mp16_3"    
}


# 2) ensure you have the subcluster object in memory as `adata_sub`
#    (whichever name you used when running the subclustering)
# If not already in RAM, load it:
adata_sub = sc.read_h5ad(r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad")

# 3) call the plotter (writes a PDF and prints the path)
pdf_path = plot_sub_umap_by_animal_mapping(
    adata_sub,               # <- your subset AnnData
    dataset_to_animal,
    datasets_col="Datasets", # or "section_id" if that's what you have
    out_dir=OUT_DIR,
    filename="sub_umap_by_animal_id.pdf",
    point_size=3.0
)
print("Saved to:", pdf_path)


#%% Plot gene expression on the subcluster UMAP by experimental group
# ==== Gene expression UMAPs, split by animal groups (PDFs) ====
# ==== Gene expression UMAPs by animal groups, with percentile color scale ====
from pathlib import PureWindowsPath, Path
from typing import Iterable, Union, Optional, List, Tuple
import numpy as np, pandas as pd, scanpy as sc, matplotlib.pyplot as plt
from scipy import sparse as sp
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def gray_bordeaux_cmap():
    # light gray → mid gray → bordeaux → very dark bordeaux
    colors = ["#E6E6E6", "#BDBDBD", "#800020", "#2B000F"]
    return LinearSegmentedColormap.from_list("gray_bordeaux", colors, N=256)


# crisp PDFs; keep points rasterized inside vector PDF
sc.settings.set_figure_params(dpi=150, dpi_save=1200, vector_friendly=True)
RASTERIZE_POINTS = True

def _canon(s):
    return str(s).strip().strip('"').strip("'").replace("\\", "/").lower()
def _name(s):  return PureWindowsPath(str(s)).name.lower()
def _stem(s):  return PureWindowsPath(str(s)).stem.lower()

def _map_animals_from_datasets(ds_series: pd.Series, mapping: dict) -> pd.Series:
    full = {_canon(k): v for k, v in mapping.items()}
    name = {_name(k):  v for k, v in mapping.items()}
    stem = {_stem(k):  v for k, v in mapping.items()}
    ds_full = ds_series.astype(str).map(_canon)
    animals = ds_full.map(full)
    miss = animals.isna()
    if miss.any(): animals.loc[miss] = ds_full[miss].map(_name).map(name)
    miss = animals.isna()
    if miss.any(): animals.loc[miss] = ds_full[miss].map(_stem).map(stem)
    return animals.fillna("UNKNOWN")

def _animal_group(animal_id: str) -> str:
    aid = (animal_id or "").lower()
    if aid.startswith("wt"): return "control_wt"
    if aid.startswith(("mp8","mt8")): return "mitopark_8w"
    if aid.startswith(("mp11","mt11")): return "mitopark_11w"
    if aid.startswith(("mp16","mt16")): return "mitopark_16w"
    return "other"

def _auto_point_size(n_obs: int) -> float:
    if n_obs <= 2_000: return 12
    if n_obs <= 10_000: return 8
    if n_obs <= 50_000: return 6
    return max(2, min(3.5, 200000 / max(1, n_obs)))

def _find_gene_case(adata, g) -> Optional[str]:
    lut = {str(v).lower(): str(v) for v in adata.var_names}
    return lut.get(str(g).lower())

def _percentile_numbers(arr, p1: int, p2: int) -> Tuple[float, float]:
    if sp.issparse(arr): arr = arr.toarray()
    x = np.asarray(arr).ravel()
    vmin, vmax = np.nanpercentile(x, [p1, p2])
    return float(vmin), float(vmax)

def plot_gene_umaps_by_animal_groups(
    adata,
    genes: Iterable[str],
    dataset_to_animal: dict,
    datasets_col: str = "Datasets",
    out_dir: Union[str, Path] = "F:/vizgen_clustering_out/subcluster_TH_DAT/umap_genes_by_group",
    file_prefix: str = "umap_expr",
    point_size: Optional[float] = None,
    percentile_scale: Tuple[str, str] = ("p1", "p99"),  # exactly like your example
    global_scale: bool = False,  # set True to use a *shared* numeric scale per gene across all groups
) -> List[Path]:
    """
    Writes up to three PDFs per gene:
      <prefix>__<gene>__control_wt.pdf
      <prefix>__<gene>__mitopark_11w.pdf
      <prefix>__<gene>__mitopark_16w.pdf

    If global_scale=False: uses Scanpy's string percentiles (e.g. vmin='p1', vmax='p99')
    per group (like your example). If True: computes numeric p1/p99 across ALL cells
    once per gene and applies the same [vmin, vmax] to each group.
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # choose dataset column
    if datasets_col not in adata.obs:
        if "section_id" in adata.obs:
            datasets_col = "section_id"
        else:
            raise KeyError(f"'{datasets_col}' (or 'section_id') not found in .obs")

    # temp group column
    animals = _map_animals_from_datasets(adata.obs[datasets_col], dataset_to_animal)
    tmp_group = "__tmp_animal_group__"
    while tmp_group in adata.obs.columns: tmp_group += "_x"
    adata.obs[tmp_group] = animals.map(_animal_group).astype("category")

    # ensure UMAP
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    s = _auto_point_size(adata.n_obs) if point_size is None else float(point_size)
    saved: List[Path] = []
    group_order = ["control_wt", "mitopark_8w", "mitopark_11w", "mitopark_16w"]

    # parse requested percentile strings once
    def _parse_p(s: str) -> int:
        return int(str(s).lower().lstrip("p"))
    p_lo = _parse_p(percentile_scale[0])
    p_hi = _parse_p(percentile_scale[1])

    for g in genes:
        gene = _find_gene_case(adata, g)
        if gene is None:
            print(f"[WARN] gene '{g}' not found; skipping.")
            continue

        # numeric global scale if requested
        vmin_num = vmax_num = None
        if global_scale:
            vmin_num, vmax_num = _percentile_numbers(adata[:, gene].X, p_lo, p_hi)

        for grp in group_order:
            mask = (adata.obs[tmp_group] == grp)
            n = int(mask.sum())
            if n == 0:
                print(f"[info] group '{grp}' has 0 cells; skipping.")
                continue

            ad_sub = adata[mask].copy()

            # choose vmin/vmax args
            if global_scale:
                vmin_arg, vmax_arg = vmin_num, vmax_num
            else:
                vmin_arg, vmax_arg = percentile_scale  # e.g. "p1","p99"

            # draw (string percentiles first; fall back to numeric if Scanpy version lacks it)
            try:
               fig = sc.pl.umap(
                        ad_sub,
                        color=gene,
                        vmin=vmin_arg, vmax=vmax_arg,   # or ("p1","p99") if you're using the percentile strings
                        s=s,
                        cmap=gray_bordeaux_cmap(),      # ← here
                        show=False,
                        return_fig=True,
                    )

            except Exception:
                # fallback: compute per-group numeric p1/p99 and retry
                if not global_scale and isinstance(vmin_arg, str) and isinstance(vmax_arg, str):
                    vmin_f, vmax_f = _percentile_numbers(ad_sub[:, gene].X, p_lo, p_hi)
                    fig = sc.pl.umap(
                        ad_sub,
                        color=gene,
                        vmin=vmin_f, vmax=vmax_f,
                        s=s,
                        cmap=gray_bordeaux_cmap(),
                        show=False,
                        return_fig=True,
                    )
                else:
                    raise

            if RASTERIZE_POINTS:
                for ax in fig.axes:
                    for coll in getattr(ax, "collections", []):
                        try: coll.set_rasterized(True)
                        except Exception: pass

            pdf_path = out_dir / f"{file_prefix}__{gene}__{grp}.pdf"
            fig.suptitle(f"{gene} • {grp} (n={n})", y=0.98, fontsize=11)
            fig.savefig(pdf_path, bbox_inches="tight")
            plt.close(fig)
            print(f"[save] {pdf_path}")
            saved.append(pdf_path)

    # cleanup temp
    del adata.obs[tmp_group]
    return saved


GENES = ["Anxa1"]  # whatever you want to inspect

#"Anxa1","Calb1", "Postn","Aldh1a1","Ndnf","Cck", "Otx2","Atf3", "Mustn1", "Tnfrsf12a","Nrp1","Epas1","Cstb","Chrm2", 
#  "Htr2a", "Grin2c", "Aldh1a7", "Atf3", "Mustn1", "Atp8b1", "Creb5", 
#"Cstb", "Cyp26b1", "Epas1", "Fkbp1b", "Chrna3", "Vip", "Egfr", "Fgf10"
#         "Hpgd", "Htr4", "Kcnj6"
#"Nr4a2", "Rmst", "Sncg", "Sprr1a", "Srpx2", "Tafa4", "Tnfrsf12a", "Zfp536","Egfr", "Fgf10", "Nrp1"

# Use your subcluster AnnData if that's the embedding you care about:
adata_sub = sc.read_h5ad(r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad")


paths = plot_gene_umaps_by_animal_groups(
    adata_sub,                  # or 'adata' if you want the full object
    GENES,
    dataset_to_animal,
    datasets_col="Datasets",    # or "section_id"
    #out_dir="F:/vizgen_clustering_out/subcluster_TH_DAT/umap_genes_by_group",
    file_prefix="sub_umap_expr",
    point_size=5.0,
    percentile_scale=("p1","p99"),  # ← matches your example
    global_scale=False               # ← per-group percentiles
)
paths = plot_gene_umaps_by_animal_groups(
    adata_sub, GENES, dataset_to_animal,
    global_scale=True, percentile_scale=("p1","p99")
)






