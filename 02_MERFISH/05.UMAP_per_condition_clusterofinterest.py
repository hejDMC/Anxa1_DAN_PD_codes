# ==== UMAPs by animal groups — CLUSTER VERSION (all-clusters + single-cluster highlight) ====
# Drop-in for your pipeline. Works with your dataset_to_animal mapping and "Datasets"/"section_id".
# Save as e.g. `umap_by_group_clusters.py` and import, or paste into your notebook.

from pathlib import Path, PureWindowsPath
from typing import Optional, Union, Iterable, List
import numpy as np, pandas as pd, scanpy as sc, matplotlib.pyplot as plt
from scipy import sparse as sp
from matplotlib.colors import LinearSegmentedColormap

# ---------- Style / shared helpers ----------
def gray_bordeaux_cmap():
    # light gray → mid gray → bordeaux → very dark bordeaux
    colors = ["#E6E6E6", "#BDBDBD", "#800020", "#2B000F"]
    return LinearSegmentedColormap.from_list("gray_bordeaux", colors, N=256)

sc.settings.set_figure_params(dpi=150, dpi_save=1200, vector_friendly=True)
RASTERIZE_POINTS = True

def _canon(s):  return str(s).strip().strip('"').strip("'").replace("\\", "/").lower()
def _name(s):   return PureWindowsPath(str(s)).name.lower()
def _stem(s):   return PureWindowsPath(str(s)).stem.lower()

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

def _match_cluster_label(values: Iterable, target) -> Optional[object]:
    """Robustly match a provided cluster label (int/str) to what's in .obs[cluster_col]."""
    if target is None:
        return None
    t_str = str(target).strip().lower()
    lut = {str(v).strip().lower(): v for v in pd.Index(values).unique()}
    if t_str in lut:
        return lut[t_str]
    try:
        t_num = int(float(t_str))
        for v in pd.Index(values).unique():
            try:
                if int(float(str(v))) == t_num:
                    return v
            except Exception:
                continue
    except Exception:
        pass
    return None

# ---------- NEW: choose the correct subclustering column automatically ----------
def _n_unique_safe(s):
    try: return int(pd.Series(s).nunique(dropna=False))
    except Exception: return 0

def _guess_cluster_col(adata, preferred: Optional[str] = None, min_levels: int = 2) -> str:
    """
    Choose a likely subclustering column from .obs.
    Priority:
      1) `preferred` if it exists and has >= min_levels unique values
      2) columns matching common subcluster patterns with >= min_levels
      3) any categorical/object/int column with >= min_levels, sorted by (#levels desc)
    """
    if preferred in adata.obs and _n_unique_safe(adata.obs[preferred]) >= min_levels:
        print(f"[cluster_col] Using preferred: '{preferred}' "
              f"({ _n_unique_safe(adata.obs[preferred]) } levels)")
        return preferred

    obs = adata.obs
    names = list(obs.columns)

    patterns_priority = [
        "subcluster", "sub_cluster", "subclust", "sub",
        "leiden_sub", "sub_leiden", "leiden_th", "leiden_dat",
        "refined", "fine", "resolution", "res", "rf", "subid", "clusters"
    ]

    def score(name: str) -> int:
        ln = name.lower()
        sc = 0
        for i, p in enumerate(patterns_priority):
            if p in ln: sc += 1000 - i
        if "leiden" in ln: sc += 50
        if "louvain" in ln: sc += 25
        return sc

    candidates = []
    for c in names:
        col = obs[c]
        if str(col.dtype) in ("category", "object") or pd.api.types.is_integer_dtype(col.dtype):
            k = _n_unique_safe(col)
            if k >= min_levels:
                candidates.append((score(c), k, c))

    if not candidates:
        raise KeyError("Could not find a cluster column with ≥2 unique levels. "
                       "Pass cluster_col= explicitly (e.g., 'sub_leiden').")

    candidates.sort(key=lambda t: (-t[0], -t[1], t[2]))

    print("[cluster_col] Candidates (score, levels, name):")
    for scv in candidates[:10]:
        print("   ", scv)

    chosen = candidates[0][2]
    print(f"[cluster_col] Auto-selected: '{chosen}' "
          f"({ _n_unique_safe(obs[chosen]) } levels)")
    return chosen

# ---------- 1) All clusters per condition ----------
def plot_cluster_umaps_by_animal_groups(
    adata,
    dataset_to_animal: dict,
    cluster_col: Optional[str] = "leiden",   # can be wrong/parent; we auto-fix
    datasets_col: str = "Datasets",
    out_dir: Union[str, Path] = "F:/vizgen_clustering_out/subcluster_TH_DAT/umap_clusters_by_group",
    file_prefix: str = "umap_clusters",
    point_size: Optional[float] = None,
    palette: Optional[Iterable[str]] = None,
    legend_loc: Optional[str] = None
) -> List[Path]:
    """
    Saves PDFs colored by clusters for each group:
      <prefix>__<group>.pdf
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # dataset column
    if datasets_col not in adata.obs:
        if "section_id" in adata.obs:
            datasets_col = "section_id"
        else:
            raise KeyError(f"'{datasets_col}' (or 'section_id') not found in .obs")

    # group assignment
    animals = _map_animals_from_datasets(adata.obs[datasets_col], dataset_to_animal)
    tmp_group = "__tmp_animal_group__"
    while tmp_group in adata.obs.columns: tmp_group += "_x"
    adata.obs[tmp_group] = animals.map(_animal_group).astype("category")

    # choose subcluster column (auto-fixes constant 'leiden')
    if (cluster_col not in adata.obs) or (_n_unique_safe(adata.obs[cluster_col]) < 2):
        cluster_col = _guess_cluster_col(adata, preferred=cluster_col, min_levels=2)
    else:
        print(f"[cluster_col] Using provided: '{cluster_col}' "
              f"({ _n_unique_safe(adata.obs[cluster_col]) } levels)")

    # UMAP
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    s = _auto_point_size(adata.n_obs) if point_size is None else float(point_size)
    group_order = ["control_wt", "mitopark_8w", "mitopark_11w", "mitopark_16w"]
    saved: List[Path] = []

    for grp in group_order:
        mask = (adata.obs[tmp_group] == grp)
        n = int(mask.sum())
        if n == 0:
            print(f"[info] group '{grp}' has 0 cells; skipping.")
            continue

        ad_sub = adata[mask].copy()

        fig = sc.pl.umap(
            ad_sub,
            color=cluster_col,
            s=s,
            palette=palette,
            legend_loc=legend_loc,
            show=False,
            return_fig=True,
        )

        if RASTERIZE_POINTS:
            for ax in fig.axes:
                for coll in getattr(ax, "collections", []):
                    try: coll.set_rasterized(False)
                    except Exception: pass

        pdf_path = out_dir / f"{file_prefix}__{grp}.pdf"
        fig.suptitle(f"{cluster_col} • {grp} (n={n})", y=0.98, fontsize=11)
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {pdf_path}")
        saved.append(pdf_path)

    del adata.obs[tmp_group]
    return saved

# ---------- 2) Single cluster highlight per condition ----------
def plot_single_cluster_highlight_umaps(
    adata,
    dataset_to_animal: dict,
    cluster_col: Optional[str] = "leiden",   # can be wrong/parent; we auto-fix
    cluster_of_interest: Union[str, int] = None,
    datasets_col: str = "Datasets",
    out_dir: Union[str, Path] = "F:/vizgen_clustering_out/subcluster_TH_DAT/umap_cluster_highlight_by_group",
    file_prefix: str = "umap_clusterHI",
    point_size: Optional[float] = None,
    fg_bg_cmap = gray_bordeaux_cmap,
) -> List[Path]:
    """
    Saves PDFs with a single cluster highlighted in bordeaux; others in gray:
      <prefix>__<clusterCol-label>__<group>.pdf
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    if datasets_col not in adata.obs:
        if "section_id" in adata.obs:
            datasets_col = "section_id"
        else:
            raise KeyError(f"'{datasets_col}' (or 'section_id') not found in .obs")

    # group assignment
    animals = _map_animals_from_datasets(adata.obs[datasets_col], dataset_to_animal)
    tmp_group = "__tmp_animal_group__"
    while tmp_group in adata.obs.columns: tmp_group += "_x"
    adata.obs[tmp_group] = animals.map(_animal_group).astype("category")

    # choose subcluster column (auto-fixes constant 'leiden')
    if (cluster_col not in adata.obs) or (_n_unique_safe(adata.obs[cluster_col]) < 2):
        cluster_col = _guess_cluster_col(adata, preferred=cluster_col, min_levels=2)
    else:
        print(f"[cluster_col] Using provided: '{cluster_col}' "
              f"({ _n_unique_safe(adata.obs[cluster_col]) } levels)")

    # label to highlight
    label = _match_cluster_label(adata.obs[cluster_col], cluster_of_interest)
    if label is None:
        uniq = pd.Index(adata.obs[cluster_col]).unique().tolist()
        raise ValueError(f"Cluster '{cluster_of_interest}' not found in '{cluster_col}'. "
                         f"Available: {uniq}")

    # UMAP
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    s = _auto_point_size(adata.n_obs) if point_size is None else float(point_size)
    group_order = ["control_wt", "mitopark_8w","mitopark_11w", "mitopark_16w"]
    saved: List[Path] = []

    tmp_bin = "__tmp_cluster_hi__"
    while tmp_bin in adata.obs.columns: tmp_bin += "_x"
    adata.obs[tmp_bin] = (adata.obs[cluster_col] == label)

    for grp in group_order:
        mask = (adata.obs[tmp_group] == grp)
        n = int(mask.sum())
        if n == 0:
            print(f"[info] group '{grp}' has 0 cells; skipping.")
            continue

        ad_sub = adata[mask].copy()

        # map True/False → 1/0 so we can apply a 2-tone cmap
        vals = ad_sub.obs[tmp_bin].astype(int)
        ad_sub.obs[tmp_bin] = vals

        fig = sc.pl.umap(
            ad_sub,
            color=tmp_bin,
            vmin=0, vmax=1,
            s=s,
            cmap=fg_bg_cmap(),
            show=False,
            return_fig=True,
        )

        if RASTERIZE_POINTS:
            for ax in fig.axes:
                for coll in getattr(ax, "collections", []):
                    try: coll.set_rasterized(True)
                    except Exception: pass

        n_fg = int(vals.sum())
        pdf_path = out_dir / f"{file_prefix}__{cluster_col}-{label}__{grp}.pdf"
        fig.suptitle(f"{cluster_col} = {label} • {grp} (n={n}; highlighted={n_fg})",
                     y=0.98, fontsize=11)
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {pdf_path}")
        saved.append(pdf_path)

    # cleanup
    del adata.obs[tmp_group]
    del adata.obs[tmp_bin]
    return saved

# ---------- Example usage ----------
if __name__ == "__main__":
    # Your mapping from the message (paste exactly as-is)
    
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
    

    # Load your subcluster AnnData
    adata_sub = sc.read_h5ad(r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad")

    # # --- A) All clusters per condition (auto-detect subcluster column if needed)
    paths_all = plot_cluster_umaps_by_animal_groups(
        adata_sub,
        dataset_to_animal=dataset_to_animal,
        cluster_col="leiden",              # if this is parent/constant, auto-picker will switch to subcluster column
        datasets_col="Datasets",           # or "section_id"
        out_dir="F:/vizgen_clustering_out/subcluster_TH_DAT/umap_clusters_by_group",
        file_prefix="sub_umap_clusters",
        point_size=5.0,
        legend_loc="right margin"
    )

    # # --- B) Highlight ONE cluster per condition (e.g., subcluster '0')
    # paths_hi = plot_single_cluster_highlight_umaps(
    #     adata_sub,
    #     dataset_to_animal=dataset_to_animal,
    #     cluster_col="leiden_sub",              # deliberately left as 'leiden'; will auto-switch if needed
    #     cluster_of_interest=1,             # accepts 0/"0"/"C0"/etc. if it matches numerically/string-wise
    #     datasets_col="Datasets",
    #     out_dir="F:/vizgen_clustering_out/subcluster_TH_DAT/umap_cluster_highlight_by_group",
    #     file_prefix="sub_umap_clusterHI",
    #     point_size=5.0
    # )

    print("\nAll-cluster PDFs:", *paths_all, sep="\n  ")
  #  print("\nHighlight PDFs:", *paths_hi, sep="\n  ")

#%% change colour of clusters in wt 



from pathlib import Path, PureWindowsPath
from typing import Optional, Union, Iterable, List
import numpy as np, pandas as pd, scanpy as sc, matplotlib.pyplot as plt
from scipy import sparse as sp
from matplotlib.colors import LinearSegmentedColormap

# ---------- Style / shared helpers ----------
def gray_bordeaux_cmap():
    # light gray → mid gray → bordeaux → very dark bordeaux
    colors = ["#E6E6E6", "#BDBDBD", "#800020", "#2B000F"]
    return LinearSegmentedColormap.from_list("gray_bordeaux", colors, N=256)

sc.settings.set_figure_params(dpi=150, dpi_save=1200, vector_friendly=True)
RASTERIZE_POINTS = True

def _canon(s):  return str(s).strip().strip('"').strip("'").replace("\\", "/").lower()
def _name(s):   return PureWindowsPath(str(s)).name.lower()
def _stem(s):   return PureWindowsPath(str(s)).stem.lower()

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

def _match_cluster_label(values: Iterable, target) -> Optional[object]:
    """Robustly match a provided cluster label (int/str) to what's in .obs[cluster_col]."""
    if target is None:
        return None
    t_str = str(target).strip().lower()
    lut = {str(v).strip().lower(): v for v in pd.Index(values).unique()}
    if t_str in lut:
        return lut[t_str]
    try:
        t_num = int(float(t_str))
        for v in pd.Index(values).unique():
            try:
                if int(float(str(v))) == t_num:
                    return v
            except Exception:
                continue
    except Exception:
        pass
    return None


def override_cluster_colors(adata, cluster_col, overrides):
    """
    overrides: dict like {'0': '#1f77b4', 5: '#d62728'}
    Only changes the specified clusters; keeps all other colors as-is.
    """
    # Ensure categorical
    adata.obs[cluster_col] = adata.obs[cluster_col].astype('category')
    cats = list(adata.obs[cluster_col].cat.categories)

    # Make sure a palette exists; if not, let scanpy create one quickly
    key = f"{cluster_col}_colors"
    if key not in adata.uns or len(adata.uns[key]) != len(cats):
        _ = sc.pl.umap(adata, color=cluster_col, show=False)  # assigns colors
        plt.close()  # don’t display

    # Build mapping cat -> current color
    current = list(adata.uns[key])
    lut = {c: current[i] for i, c in enumerate(cats)}

    # Apply overrides (robust to int/str labels)
    for k, new_col in overrides.items():
        # match label as string or numeric
        targets = [lab for lab in cats
                   if str(lab) == str(k) or (str(lab).isdigit() and str(k).isdigit() and int(lab) == int(k))]
        for t in targets:
            lut[t] = new_col

    # Write back in category order
    adata.uns[key] = [lut[c] for c in cats]


# ---------- NEW: choose the correct subclustering column automatically ----------
def _n_unique_safe(s):
    try: return int(pd.Series(s).nunique(dropna=False))
    except Exception: return 0

def _guess_cluster_col(adata, preferred: Optional[str] = None, min_levels: int = 2) -> str:
    """
    Choose a likely subclustering column from .obs.
    Priority:
      1) `preferred` if it exists and has >= min_levels unique values
      2) columns matching common subcluster patterns with >= min_levels
      3) any categorical/object/int column with >= min_levels, sorted by (#levels desc)
    """
    if preferred in adata.obs and _n_unique_safe(adata.obs[preferred]) >= min_levels:
        print(f"[cluster_col] Using preferred: '{preferred}' "
              f"({ _n_unique_safe(adata.obs[preferred]) } levels)")
        return preferred

    obs = adata.obs
    names = list(obs.columns)

    patterns_priority = [
        "subcluster", "sub_cluster", "subclust", "sub",
        "leiden_sub", "sub_leiden", "leiden_th", "leiden_dat",
        "refined", "fine", "resolution", "res", "rf", "subid", "clusters"
    ]

    def score(name: str) -> int:
        ln = name.lower()
        sc = 0
        for i, p in enumerate(patterns_priority):
            if p in ln: sc += 1000 - i
        if "leiden" in ln: sc += 50
        if "louvain" in ln: sc += 25
        return sc

    candidates = []
    for c in names:
        col = obs[c]
        if str(col.dtype) in ("category", "object") or pd.api.types.is_integer_dtype(col.dtype):
            k = _n_unique_safe(col)
            if k >= min_levels:
                candidates.append((score(c), k, c))

    if not candidates:
        raise KeyError("Could not find a cluster column with ≥2 unique levels. "
                       "Pass cluster_col= explicitly (e.g., 'sub_leiden').")

    candidates.sort(key=lambda t: (-t[0], -t[1], t[2]))

    print("[cluster_col] Candidates (score, levels, name):")
    for scv in candidates[:10]:
        print("   ", scv)

    chosen = candidates[0][2]
    print(f"[cluster_col] Auto-selected: '{chosen}' "
          f"({ _n_unique_safe(obs[chosen]) } levels)")
    return chosen



cluster_col = "leiden_sub"  # or the one you actually plot with (your auto-guesser may change it)
override_cluster_colors(
    adata_sub,
    cluster_col,
    overrides={'0': '#993a52', '8': '#8e91a7'}  # pick your two clusters and colors
)

# then run your existing plot function(s)
paths_all = plot_cluster_umaps_by_animal_groups(
    adata_sub,
    dataset_to_animal=dataset_to_animal,
    cluster_col=cluster_col,
    datasets_col="Datasets",
    out_dir="F:/vizgen_clustering_out/subcluster_TH_DAT/umap_clusters_by_group",
    file_prefix="sub_umap_clusters",
    point_size=5.0,
    legend_loc="right margin"
)




# ---------- 1) All clusters per condition ----------
def plot_cluster_umaps_by_animal_groups(
    adata,
    dataset_to_animal: dict,
    cluster_col: Optional[str] = "leiden",   # can be wrong/parent; we auto-fix
    datasets_col: str = "Datasets",
    out_dir: Union[str, Path] = "F:/vizgen_clustering_out/subcluster_TH_DAT/umap_clusters_by_group",
    file_prefix: str = "umap_clusters",
    point_size: Optional[float] = None,
    palette: Optional[Iterable[str]] = None,
    legend_loc: Optional[str] = None
) -> List[Path]:
    """
    Saves PDFs colored by clusters for each group:
      <prefix>__<group>.pdf
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # dataset column
    if datasets_col not in adata.obs:
        if "section_id" in adata.obs:
            datasets_col = "section_id"
        else:
            raise KeyError(f"'{datasets_col}' (or 'section_id') not found in .obs")

    # group assignment
    animals = _map_animals_from_datasets(adata.obs[datasets_col], dataset_to_animal)
    tmp_group = "__tmp_animal_group__"
    while tmp_group in adata.obs.columns: tmp_group += "_x"
    adata.obs[tmp_group] = animals.map(_animal_group).astype("category")

    # choose subcluster column (auto-fixes constant 'leiden')
    if (cluster_col not in adata.obs) or (_n_unique_safe(adata.obs[cluster_col]) < 2):
        cluster_col = _guess_cluster_col(adata, preferred=cluster_col, min_levels=2)
    else:
        print(f"[cluster_col] Using provided: '{cluster_col}' "
              f"({ _n_unique_safe(adata.obs[cluster_col]) } levels)")

    # UMAP
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    s = _auto_point_size(adata.n_obs) if point_size is None else float(point_size)
    group_order = ["control_wt", "mitopark_8w", "mitopark_11w", "mitopark_16w"]
    saved: List[Path] = []

    for grp in group_order:
        mask = (adata.obs[tmp_group] == grp)
        n = int(mask.sum())
        if n == 0:
            print(f"[info] group '{grp}' has 0 cells; skipping.")
            continue

        ad_sub = adata[mask].copy()

        fig = sc.pl.umap(
            ad_sub,
            color=cluster_col,
            s=s,
            palette=palette,
            legend_loc=legend_loc,
            show=False,
            return_fig=True,
        )

        if RASTERIZE_POINTS:
            for ax in fig.axes:
                for coll in getattr(ax, "collections", []):
                    try: coll.set_rasterized(False)
                    except Exception: pass

        pdf_path = out_dir / f"{file_prefix}__{grp}.pdf"
        fig.suptitle(f"{cluster_col} • {grp} (n={n})", y=0.98, fontsize=11)
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {pdf_path}")
        saved.append(pdf_path)

    del adata.obs[tmp_group]
    return saved

# ---------- 2) Single cluster highlight per condition ----------
def plot_single_cluster_highlight_umaps(
    adata,
    dataset_to_animal: dict,
    cluster_col: Optional[str] = "leiden",   # can be wrong/parent; we auto-fix
    cluster_of_interest: Union[str, int] = None,
    datasets_col: str = "Datasets",
    out_dir: Union[str, Path] = "F:/vizgen_clustering_out/subcluster_TH_DAT/umap_cluster_highlight_by_group",
    file_prefix: str = "umap_clusterHI",
    point_size: Optional[float] = None,
    fg_bg_cmap = gray_bordeaux_cmap,
) -> List[Path]:
    """
    Saves PDFs with a single cluster highlighted in bordeaux; others in gray:
      <prefix>__<clusterCol-label>__<group>.pdf
    """
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    if datasets_col not in adata.obs:
        if "section_id" in adata.obs:
            datasets_col = "section_id"
        else:
            raise KeyError(f"'{datasets_col}' (or 'section_id') not found in .obs")

    # group assignment
    animals = _map_animals_from_datasets(adata.obs[datasets_col], dataset_to_animal)
    tmp_group = "__tmp_animal_group__"
    while tmp_group in adata.obs.columns: tmp_group += "_x"
    adata.obs[tmp_group] = animals.map(_animal_group).astype("category")

    # choose subcluster column (auto-fixes constant 'leiden')
    if (cluster_col not in adata.obs) or (_n_unique_safe(adata.obs[cluster_col]) < 2):
        cluster_col = _guess_cluster_col(adata, preferred=cluster_col, min_levels=2)
    else:
        print(f"[cluster_col] Using provided: '{cluster_col}' "
              f"({ _n_unique_safe(adata.obs[cluster_col]) } levels)")

    # label to highlight
    label = _match_cluster_label(adata.obs[cluster_col], cluster_of_interest)
    if label is None:
        uniq = pd.Index(adata.obs[cluster_col]).unique().tolist()
        raise ValueError(f"Cluster '{cluster_of_interest}' not found in '{cluster_col}'. "
                         f"Available: {uniq}")

    # UMAP
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    s = _auto_point_size(adata.n_obs) if point_size is None else float(point_size)
    group_order = ["control_wt", "mitopark_8w","mitopark_11w", "mitopark_16w"]
    saved: List[Path] = []

    tmp_bin = "__tmp_cluster_hi__"
    while tmp_bin in adata.obs.columns: tmp_bin += "_x"
    adata.obs[tmp_bin] = (adata.obs[cluster_col] == label)

    for grp in group_order:
        mask = (adata.obs[tmp_group] == grp)
        n = int(mask.sum())
        if n == 0:
            print(f"[info] group '{grp}' has 0 cells; skipping.")
            continue

        ad_sub = adata[mask].copy()

        # map True/False → 1/0 so we can apply a 2-tone cmap
        vals = ad_sub.obs[tmp_bin].astype(int)
        ad_sub.obs[tmp_bin] = vals

        fig = sc.pl.umap(
            ad_sub,
            color=tmp_bin,
            vmin=0, vmax=1,
            s=s,
            cmap=fg_bg_cmap(),
            show=False,
            return_fig=True,
        )

        if RASTERIZE_POINTS:
            for ax in fig.axes:
                for coll in getattr(ax, "collections", []):
                    try: coll.set_rasterized(True)
                    except Exception: pass

        n_fg = int(vals.sum())
        pdf_path = out_dir / f"{file_prefix}__{cluster_col}-{label}__{grp}.pdf"
        fig.suptitle(f"{cluster_col} = {label} • {grp} (n={n}; highlighted={n_fg})",
                     y=0.98, fontsize=11)
        fig.savefig(pdf_path, bbox_inches="tight")
        plt.close(fig)
        print(f"[save] {pdf_path}")
        saved.append(pdf_path)

    # cleanup
    del adata.obs[tmp_group]
    del adata.obs[tmp_bin]
    return saved

# ---------- Example usage ----------
if __name__ == "__main__":
    # Your mapping from the message (paste exactly as-is)
    
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
    

    # Load your subcluster AnnData
    adata_sub = sc.read_h5ad(r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad")

    # # --- A) All clusters per condition (auto-detect subcluster column if needed)
    paths_all = plot_cluster_umaps_by_animal_groups(
        adata_sub,
        dataset_to_animal=dataset_to_animal,
        cluster_col="leiden",              # if this is parent/constant, auto-picker will switch to subcluster column
        datasets_col="Datasets",           # or "section_id"
        out_dir="F:/vizgen_clustering_out/subcluster_TH_DAT/umap_clusters_by_group",
        file_prefix="sub_umap_clusters",
        point_size=5.0,
        legend_loc="right margin"
    )


    print("\nAll-cluster PDFs:", *paths_all, sep="\n  ")
  #  print("\nHighlight PDFs:", *paths_hi, sep="\n  ")






















