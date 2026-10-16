import scanpy as sc
import numpy as np
from pathlib import Path

def dotplot_selected_genes(
    ad,
    genes,                      # list of gene symbols in the exact order you want
    cluster_key,                # e.g. "leiden_sub_anxa1pos" or "leiden_sub"
    layer="log1p",              # or "logexpr" or None to use .X
    cmap="Reds",
    standard_scale="var",       # None for absolute colors
    swap_axes=False,            # True → genes on y-axis
    out_pdf=r"dotplot_selected_genes.pdf",
):
    # case-insensitive match to ad.var_names, preserve requested order
    look = {g.lower(): g for g in ad.var_names}
    found, missing = [], []
    for g in genes:
        gg = look.get(str(g).lower())
        (found.append(gg) if gg is not None else missing.append(g))
    if missing:
        print(f"[warn] not found: {', '.join(map(str, missing))}")
    if not found:
        raise ValueError("None of the requested genes were found in var_names.")

    dp = sc.pl.dotplot(
        ad,
        var_names=found,               # exact genes you picked, ordered
        groupby=cluster_key,
        layer=layer if layer is not None else None,
        use_raw=False,
        standard_scale=standard_scale,
        dot_min=0.0, dot_max=1.0,
        cmap=cmap,                     # Scanpy 1.9.3 uses 'cmap'
        swap_axes=swap_axes,
        return_fig=True, show=False
    )
    # keep the same colormap after styling (1.9.3 quirk)
    dp = dp.style(cmap=cmap, dot_edge_color=None, dot_edge_lw=0.0)
    Path(out_pdf).parent.mkdir(parents=True, exist_ok=True)
    dp.savefig(out_pdf)
    print(f"[saved] {out_pdf}")
    return out_pdf


SUBCLUSTER_H5AD = r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad"
ad = sc.read_h5ad(SUBCLUSTER_H5AD)
CLUST_KEY = "leiden_sub_anxa1pos" if "leiden_sub_anxa1pos" in ad.obs else "leiden_sub"

MARKERS = ["Cyp26b1","Tacr3",'Ndnf',"Igf1", #0
           'Aldh1a7', 'Anxa1','Aldh1a1',  "Srpx2", #1
           "Ebf2", "Npnt", #2
           "Atp2a3", 
           "Npnt", "Fst","Col11a1" ,"Col14a1", "Otx2", "Ano2",  "Rxfp1", #3
           'Sprr1a', 'Grp', "Fst","Lpl",  #4
           "Gipr","Lix1", "Gfra1", #5
           'Atf3', "Mmp12", "Vill",'Tnfrsf12a', 'Mustn1', #6
           "Calb1", #7
           "Pde11a","Pcsk6", "Dkk2" ,#8
           "Lama3",'Igf1',  "Grin2c", "Kcns3",  #9
           "Il1rapl2","Sema5b", 'Unc5d', #10
           "Wnt7b",'Postn','Npy1r',"Rxfp1", #11
           "Pdgfd","Vcan","Lix1","Kcnj6", "Grik1",'Col11a1', "Cpne2",  #12
           "Ebf1","Cck", "Slc17a6",  #13
           "Klhl14",  #14
           "Col24a1", "Fkbp1b", "Pld5", #15
           "Lpl",   #116

]  # <- YOUR list & order

dotplot_selected_genes(
    ad,
    genes=MARKERS,
    cluster_key=CLUST_KEY,
    layer="log1p" if "log1p" in ad.layers else ("logexpr" if "logexpr" in ad.layers else None),
    cmap="Reds",
    standard_scale="var",  # set None if you want absolute colors
    swap_axes=False,
    out_pdf=r"F:/vizgen_clustering_out/subcluster_TH_DAT/dotplot_selected_genes.pdf"
)

ad.obs[CLUST_KEY] = ad.obs[CLUST_KEY].cat.reorder_categories(["0","1","2","3","4","5","6","7","8","9","10"])


#%%
import scanpy as sc
import numpy as np
import pandas as pd
from pathlib import Path

def dotplot_selected_genes_grouped(
    ad,
    cluster_key,                        # e.g. "leiden_sub_anxa1pos" or "leiden_sub"
    genes_by_cluster,                   # dict: {"0":[...], "1":[...], ...} OR single list to repeat
    layer="log1p",                      # or "logexpr" or None -> use .X
    cmap="Reds",
    standard_scale="var",               # None for absolute colors
    out_pdf=r"dotplot_selected_genes_grouped.pdf",
):
    # 1) cluster order on y-axis
    if hasattr(ad.obs[cluster_key].dtype, "categories"):
        clusters = list(ad.obs[cluster_key].cat.categories.astype(str))
    else:
        clusters = sorted(ad.obs[cluster_key].astype(str).unique().tolist())
        ad.obs[cluster_key] = pd.Categorical(ad.obs[cluster_key].astype(str),
                                             categories=clusters, ordered=True)

    # 2) normalize input: if user passed a single list, repeat for each cluster label
    if isinstance(genes_by_cluster, (list, tuple)):
        genes_by_cluster = {c: list(genes_by_cluster) for c in clusters}

    # 3) case-insensitive matching of gene names, preserve given order
    look = {g.lower(): g for g in ad.var_names.astype(str)}
    var_groups, group_labels = [], []
    for c in clusters:
        wanted = genes_by_cluster.get(c, [])
        found, missing = [], []
        for g in wanted:
            gg = look.get(str(g).lower())
            (found.append(gg) if gg is not None else missing.append(g))
        if missing:
            print(f"[warn] cluster {c}: not found -> {', '.join(map(str, missing))}")
        if len(found) == 0:
            continue
        var_groups.append(found)
        group_labels.append(c)

    # 4) flatten + group positions (Scanpy 1.9.3 expects (start,end) tuples)
    flat_genes = [g for grp in var_groups for g in grp]
    lengths = [len(grp) for grp in var_groups]
    if sum(lengths) == 0:
        raise ValueError("No genes matched var_names.")
    starts = np.cumsum([0] + lengths[:-1])
    ends   = np.cumsum(lengths) - 1
    var_group_positions = list(zip(starts.tolist(), ends.tolist()))

    # 5) make the dotplot
    dp = sc.pl.dotplot(
        ad,
        var_names=flat_genes,
        var_group_labels=group_labels,          # top brackets: cluster labels
        var_group_positions=var_group_positions,
        groupby=cluster_key,                    # y-axis: clusters
        layer=layer if layer is not None else None,
        use_raw=False,
        standard_scale=standard_scale,
        dot_min=0.0, dot_max=1.0,
        cmap=cmap,                              # Scanpy 1.9.3 kw
        swap_axes=False,                        # keep genes along bottom
        return_fig=True, show=False
    )
    # keep the same cmap after styling (1.9.3 quirk)
    dp = dp.style(cmap=cmap, dot_edge_color=None, dot_edge_lw=0.0)
    Path(out_pdf).parent.mkdir(parents=True, exist_ok=True)
    dp.savefig(out_pdf)
    print(f"[saved] {out_pdf}")
    return out_pdf

SUBCLUSTER_H5AD = r"F:/vizgen_clustering_out/subcluster_TH_DAT/subset_Th_Slc6a3.h5ad"
ad = sc.read_h5ad(SUBCLUSTER_H5AD)
CLUST_KEY = "leiden_sub_anxa1pos" if "leiden_sub_anxa1pos" in ad.obs else "leiden_sub"

# fix a consistent cluster order (keep the ones that exist)
desired = [str(i) for i in range(0, 11)]
present = ad.obs[CLUST_KEY].astype(str).unique().tolist()
order = [c for c in desired if c in present]
ad.obs[CLUST_KEY] = pd.Categorical(ad.obs[CLUST_KEY].astype(str), categories=order, ordered=True)

# your exact genes & order
MARKERS = ["Aldh1a1","Anxa1","Ndnf","Postn","Vill","Calb1","Cck","Otx2",
           "Atf3","Mustn1","Tnfrsf12a","Cyp26b1","Gad2","Grin2c","Htr2a",
           "Sema5b","Slc17a6","Slc32a1"]

LAYER = "log1p" if "log1p" in ad.layers else ("logexpr" if "logexpr" in ad.layers else None)

# same markers under each top cluster bracket:
dotplot_selected_genes_grouped(
    ad,
    cluster_key=CLUST_KEY,
    genes_by_cluster=MARKERS,   # a single list → repeated for each cluster group
    layer=LAYER,
    cmap="Reds",
    standard_scale="var",       # set None for absolute colors
    out_pdf=r"F:/vizgen_clustering_out/subcluster_TH_DAT/dotplot_selected_genes__grouped_by_cluster.pdf"
)

GENES_PER_CLUSTER = {
    "0": ['Aldh1a7', 'Aldh1a1', 'Anxa1', 'Ndnf'],
    "1" : [ 'Aldh1a1', 'Anxa1','Cyp26b1', 'Ndnf','Postn', 'Vill'],
    "2" :['Gfra1', 'Atp2a3', 'Pld5', 'Col25a1'],
    "3" : ["Sema5b",'Calb1', "Npy1r", "Postn", 'Npy1r'],
    "4" : ['Grp', 'Tacr3', 'Aldh1a1', "Grp", "Lpl", 'Anxa1', 'Col11a1', "Neurod6", 'Calb1'],
    "5":['Calb1', 'Col14a1','Col11a1', "Fgf10", "Dkk2" ],
    "6": ['Vill', 'Ndnf', 'Igf1', 'Lmcd1',  "Grin2c"],
    "7" : ['Atf3', 'Vill', 'Tnfrsf12a', 'Mustn1', 'Anxa1', 'Cstb', "Aldh1a1"],
    "8": ['Calb1', 'Cyp26b1', 'Col24a1', "Vip", 'Unc5d'],
    "9" : ['Vcan', 'Pld5', 'Ebf2', 'Npnt'],
    "10": ['Vcan', 'Atf3', 'Htr4', ]
    # fill others; any clusters you omit will be skipped in the top brackets
}
dotplot_selected_genes_grouped(
    ad, CLUST_KEY, GENES_PER_CLUSTER, layer=LAYER, cmap="Reds",
    out_pdf=r".../dotplot_selected_genes__per_cluster_custom.pdf"
)


#'Cyp26b1', 'Ndnf','Postn', 'Vill', 'Gfra1', 'Atp2a3', 'Pld5', 'Col25a1',"Sema5b",'Calb1', "Postn", 'Npy1r'
#'Grp', 'Tacr3', 'Aldh1a1', "Lpl", 'Anxa1', 'Col11a1', "Neurod6", 'Col14a1', "Fgf10", "Dkk2",'Vill', 'Ndnf', 'Igf1', 'Lmcd1',  "Grin2c"
#'Atf3', 'Vill', 'Tnfrsf12a', 'Mustn1', 'Anxa1', 'Cstb', "Aldh1a1", 'Calb1', 'Cyp26b1', 'Col24a1', "Vip", 'Unc5d','Vcan', 'Pld5', 'Ebf2', 'Npnt'
#'Htr4'
