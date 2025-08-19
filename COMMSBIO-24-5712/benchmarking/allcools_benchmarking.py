Python 3.10.0 (v3.10.0:b494f5935c, Oct  4 2021, 14:59:19) [Clang 12.0.5 (clang-1205.0.22.11)] on darwin
Type "help", "copyright", "credits" or "license()" for more information.

########################################################################
### Clustering with 100kb windows ###

# Benchmarking sciMETv2_3842F data for amethyst paper
# 03/25/25
# Following tutorial https://lhqing.github.io/ALLCools/cell_level/basic/mch_mcg_100k_basic.html

# cd /secret/path/amethyst/manuscript_analysis_v3/benchmark/allcools/sciMETv2_3842F_process
# wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import resource
import time

from memory_profiler import memory_usage
from ALLCools.mcds import MCDS
from ALLCools.clustering import tsne, significant_pc_test, log_scale
from ALLCools.plot import *

# Define filtering metrics
# Data already filtered; not adding additional parameters for benchmarking purposes

metadata_path = '/secret/path/amethyst/manuscript_analysis_v3/benchmark/allcools/sciMETv2_3842F_process/brain_metadata.txt'
mcds_path = '/secret/path/amethyst/manuscript_analysis_v3/benchmark/allcools/sciMETv2_3842F_process/sciMETv2_3842F.mcds'

# Dimension name used to do clustering
# This corresponding to AnnData .obs and .var
obs_dim = 'cell'  # observation
var_dim = 'chrom100k'  # feature
mch_col_name = 'mch_pct'
mch_col_name = 'mcg_pct'

# feature cov cutoffs
min_cov = 5

# Regions to remove during the clustering analysis
# change this to the path to ENCODE blacklist.
# The ENCODE blacklist can be downloaded from https://github.com/Boyle-Lab/Blacklist/
black_list_path = '/secret/path/amethyst/manuscript_analysis_v3/benchmark/allcools/sciMETv2_3842F_process/hg38-blacklist.v2.bed.gz'
black_list_fraction = 0.2
exclude_chromosome = ['chrM', 'chrY']

# load to memory or not
load = True

# HVF
mch_pattern = 'CH'
mcg_pattern = 'CG'
n_top_feature = 20000

# PC cutoff
pc_cutoff = 0.1

# KNN
knn = -1  # -1 means auto determine

# Leiden
resolution = 1

metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
print(f'Metadata of {metadata.shape[0]} cells') # 1346 cells
metadata.head()

# dataset already filtered
judge = (metadata["mch_pct"] < 12) & \
        (metadata["cov"] > 6000000) & \
        (metadata["cov"] < 120000000)

metadata = metadata[judge].copy()

# Load data
mcds = MCDS.open(
    mcds_path, 
    obs_dim='cell', 
    var_dim='chrom100k',
    use_obs=metadata.index  # MCDS contains all cells, this will select cells that passed filtering 
)
print(mcds)
total_feature = mcds.get_index(var_dim).size
mcds

mcds.add_cell_metadata(metadata)
mcds.add_feature_cov_mean(var_dim=var_dim)

# filter by coverage - based on the distribution above
mcds = mcds.filter_feature_by_cov_mean(
    min_cov=min_cov
)

# remove blacklist regions
mcds = mcds.remove_black_list_region(
    black_list_path=black_list_path,
    f=black_list_fraction  # Features having overlap > f with any black list region will be removed.
)

# remove chromosomes
mcds = mcds.remove_chromosome(exclude_chromosome)

mcds.add_mc_frac(
normalize_per_cell=True,  # after calculating mC frac, per cell normalize the matrix
    clip_norm_value=10  # clip outlier values above 10 to 10
)

# load only the mC fraction matrix into memory so following steps is faster
# Only load into memory when your memory size is enough to handle your dataset
if load and (mcds.get_index(obs_dim).size < 20000):
    mcds[f'{var_dim}_da_frac'].load()

mch_hvf = mcds.calculate_hvf_svr(var_dim=var_dim,
                                 mc_type=mch_pattern,
                                 n_top_feature=n_top_feature,
                                 plot=True)

# calculate time/memory for mCG highly variable feature selection
hvf_start_time = time.time()
mcg_hvf = mcds.calculate_hvf_svr(var_dim=var_dim,
                                 mc_type=mcg_pattern,
                                 n_top_feature=n_top_feature,
                                 plot=True)
mem_kb, result = memory_usage(
    (mcds.calculate_hvf_svr, [], {
        'var_dim': var_dim,
        'mc_type': mcg_pattern,
        'n_top_feature': n_top_feature,
        'plot': True
    }),
    retval=True,
    max_usage=True
) # Peak memory: 1887608 KB
hvf_end_time = time.time()
hvf_elapsed_minutes = (hvf_end_time - hvf_start_time) / 60 # 0.6311713218688965

mch_adata = mcds.get_adata(mc_type=mch_pattern,
                           var_dim=var_dim,
                           select_hvf=True)

mcg_adata = mcds.get_adata(mc_type=mcg_pattern,
                           var_dim=var_dim,
                           select_hvf=True)

log_scale(mch_adata)
log_scale(mcg_adata)

sc.tl.pca(mch_adata)
ch_n_components = significant_pc_test(mch_adata)


fig, axes = plot_decomp_scatters(mch_adata,	
                                 n_components=ch_n_components,
                                 hue=mch_col_name,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=3,
                                 ncols=5)

# test CG pca time and memory stats
pca_start_time = time.time()
sc.tl.pca(mcg_adata)
cg_n_components = significant_pc_test(mcg_adata)
pca_end_time = time.time()
pca_elapsed_minutes = (pca_end_time - pca_start_time) / 60 # 0.1314203143119812

def run_pca_and_test():
    sc.tl.pca(mcg_adata)
    return significant_pc_test(mcg_adata)

    
mem_kb, cg_n_components = memory_usage(run_pca_and_test, retval=True, max_usage=True) # 2442744 KB


fig, axes = plot_decomp_scatters(mcg_adata,
                                 n_components=cg_n_components,
                                 hue=mcg_col_name,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=3,
                                 ncols=5)

ch_pcs = mch_adata.obsm['X_pca'][:, :ch_n_components]
cg_pcs = mcg_adata.obsm['X_pca'][:, :cg_n_components]

# scale the PCs so CH and CG PCs has the same total var
cg_pcs = cg_pcs / cg_pcs.std()
ch_pcs = ch_pcs / ch_pcs.std()

# total_pcs
total_pcs = np.hstack([ch_pcs, cg_pcs])

# make a copy of adata, add new pcs
# this is suboptimal, will change this when adata can combine layer and X in the future
adata = mch_adata.copy()
adata.obsm['X_pca'] = total_pcs
adata.obs['type'] = metadata['type']
del adata.uns['pca']
del adata.varm['PCs']

# CLUSTERING
# Calculate nearest neighbors

if knn == -1:
    knn = max(15, int(np.log2(adata.shape[0])*2))

sc.pp.neighbors(adata, n_neighbors=knn)
sc.tl.leiden(adata, resolution=resolution)

tsne_start_time = time.time()
tsne(adata,
     obsm='X_pca',
     metric='euclidean',
     exaggeration=-1,  # auto determined
     perplexity=30,
     n_jobs=-1)

tsne_end_time = time.time()
tsne_elapsed_minutes = (tsne_end_time - tsne_start_time) / 60 # 0.4081445058186849


fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=adata,
                        ax=ax,
                        coord_base='tsne',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)

sc.tl.umap(adata)

fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=adata,
                        ax=ax,
                        coord_base='umap',
                        hue='type',
                        text_anno='type',
                        show_legend=True)
# fig.savefig("sciMETv2_3842F_umap_plot.pdf", format="pdf", bbox_inches="tight")

# SAVE RESULTS
# adata.write_h5ad('sciMETv2_3842F.chrom100k-clustering.h5ad')

########################################################################
# adapted from https://github.com/lhqing/ALLCools/blob/master/ALLCools/dmr/call_dmr.py because of source code errors

import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import xarray as xr
import time

from memory_profiler import memory_usage
from ALLCools.mcds import RegionDS
from ALLCools.mcds.utilities import update_dataset_config, write_ordered_chunks

p_value_cutoff=0.001
frac_delta_cutoff=0.2
max_dist=250
residual_quantile=0.6
corr_cutoff=0.3
dms_ratio=0.8

output_dir = "/secret/path/amethyst/manuscript_analysis_v3/benchmark/allcools/dmr/sciMETv2_3842F_pseudobulk_allc/sciMETv2_3842F_exc_vs_inh_dms"
dms_dir = "/secret/path/amethyst/manuscript_analysis_v3/benchmark/allcools/dmr/sciMETv2_3842F_pseudobulk_allc/sciMETv2_3842F_exc_vs_inh_dms/dms"
dmr_dir = "/secret/path/amethyst/manuscript_analysis_v3/benchmark/allcools/dmr/sciMETv2_3842F_pseudobulk_allc/sciMETv2_3842F_exc_vs_inh_dms/dmr"

ds = RegionDS.open(dms_dir, region_dim="dms", chunks={"dms": 1000000})

# open DMS dataset and select single chromosome, significant, large delta DMS
p_value_judge = ds.coords["dms_p-values"] < p_value_cutoff
frac_delta = ds["dms_da_frac"].max(dim="sample") - ds["dms_da_frac"].min(dim="sample")
frac_delta_judge = frac_delta > frac_delta_cutoff
# chrom_judge = ds["dms_chrom"] == chrom
# has to use load(scheduler='sync') when this function is called from multiprocess
final_judge = (p_value_judge & frac_delta_judge).to_pandas()

use_dms = final_judge[final_judge].index
ds = ds.sel(dms=use_dms)

# step 1: combine raw DMR windows based on distance
dms_dist = ds["dms_pos"][1:].values - ds["dms_pos"][:-1].values
dist_judge = dms_dist < max_dist
cur_dmr = 0
dmr_ids = [0]
for i in dist_judge:
    if not i:
        cur_dmr += 1
    dmr_ids.append(cur_dmr)
    
dmr_ids = pd.Series(dmr_ids, index=ds.get_index("dms"))

# step 2: determine correlation between adjacent DMS
data = ds["dms_da_frac"].transpose("dms", "sample").to_pandas()
a = data.iloc[:-1, :].reset_index(drop=True)
b = data.iloc[1:, :].reset_index(drop=True)
# index of corr means the correlation of that DMS with the previous DMS
# regardless of the genome distance
# fill na value with 1, tend to not to split DMR due to nan
corrs = a.fillna(1).corrwith(b.fillna(1), axis=1, method="pearson")
corrs.index = data.iloc[1:, :].index.copy()

# step 3: recalculate DMR windows based on both distance and correlation
raw_dmr_table = pd.DataFrame({"dmr": dmr_ids, "corr": corrs})
dmr_dict = {}
cur_dmr_id = 0
for _, sub_df in raw_dmr_table.groupby("dmr"):
    dmr_dict[sub_df.index[0]] = cur_dmr_id
    if sub_df.shape[0] > 1:
        for dms_id, corr in sub_df["corr"][1:].items():
            if corr > corr_cutoff:
                dmr_dict[dms_id] = cur_dmr_id
            else:
                cur_dmr_id += 1
                dmr_dict[dms_id] = cur_dmr_id
    cur_dmr_id += 1

dmr_ids = pd.Series(dmr_dict)
dmr_ids.index.name = "dms"
ds = ds.assign_coords({"dmr": ("dms", dmr_ids.reindex(ds.coords["dms"]))})

# step 4: determine sample hypo or hyper in each DMS and DMR based on residual
if residual_quantile < 0.5:
    residual_quantile = 1 - residual_quantile

residual_low_cutoff, residual_high_cutoff = np.nanquantile(
    ds["dms_residual"], [1 - residual_quantile, residual_quantile]
)

lower_residual = ds["dms_residual"] < residual_low_cutoff
higher_residual = ds["dms_residual"] > residual_high_cutoff
# for each sample in each DMS
# -1 means hypo methylation, 1 means hyper methylation, 0 mean no sig change
dms_states = lower_residual.astype(int) * -1 + higher_residual.astype(int)
# dmr state judge
# use mean to calculate the overall DMR state from DMS
dmr_states = dms_states.groupby("dmr").mean()

# then mask inconsistent DMR states with the dms_ratio cutoff
dmr_states = xr.where(np.abs(dmr_states) < dms_ratio, 0, dmr_states)
# then "round" the dmr_states, so it only contains -1, 0, 1
dmr_states = xr.where(dmr_states > 0, 1, dmr_states)
dmr_states = xr.where(dmr_states < 0, -1, dmr_states)

# step 5: prepare dmr counts and fractions
dmr_da = ds["dms_da"].groupby("dmr").sum()
dmr_frac = dmr_da.sel(count_type="mc") / dmr_da.sel(count_type="cov")
dmr_pval = ds["dms_p-values"].groupby(ds["dmr"]).mean()

dmr_ds = xr.Dataset(
    {
        "dmr_da": dmr_da.astype(np.uint32),
        "dmr_state": dmr_states.astype(np.int8),
        "dmr_da_frac": dmr_frac.astype(np.float32),
        "dmr_pval": dmr_pval.astype(np.float64)
    }
)

# add n dms counts
n_dms = dmr_ids.value_counts().sort_index()
n_dms.index.name = "dmr"
dmr_ds.coords["dmr_ndms"] = n_dms

# add genome coords
#dmr_ds.coords["dmr"] = chrom + "-" + dmr_ds["dmr"].astype(str)
dmr_ds.coords["dmr_chrom"] = ds["dms_chrom"].groupby(ds["dmr"]).first()
dmr_ds.coords["dmr_start"] = ds["dms_pos"].groupby(ds["dmr"]).min().to_pandas() - 1
dmr_ds.coords["dmr_end"] = ds["dms_pos"].groupby(ds["dmr"]).max().to_pandas() + 1

# save
dmr_ds.to_zarr(dmr_dir, mode="w")

### CONVERT TO DATA FRAME - not ALLCools code ###

# Extract only the coordinate variables first (dmr-level metadata)
df = pd.DataFrame({
    "dmr": dmr_ds["dmr"].values,  # Extract as NumPy array (1D)
    "dmr_chrom": dmr_ds["dmr_chrom"].values,  # Convert Dask array to NumPy
    "dmr_start": dmr_ds["dmr_start"].values,
    "dmr_end": dmr_ds["dmr_end"].values,
    "dmr_pval": dmr_ds["dmr_pval"].values
})

# Convert multi-sample variables while keeping one row per DMR
# Aggregating multi-sample values into separate columns
dmr_state_df = dmr_ds.dmr_state.to_pandas().T  # Transpose so rows are DMRs
dmr_da_frac_df = dmr_ds.dmr_da_frac.to_pandas().T

# Renaming columns to include sample info
dmr_state_df.columns = [f"dmr_state_{s}" for s in dmr_ds.sample.values]
dmr_da_frac_df.columns = [f"dmr_da_frac_{s}" for s in dmr_ds.sample.values]

# Merge into the final DataFrame
df = df.merge(dmr_state_df, left_on="dmr", right_index=True)
df = df.merge(dmr_da_frac_df, left_on="dmr", right_index=True)

df_filtered = df[
    ~((df["dmr_state_sciMETv2_3842F_inhibitory"] == 0) & 
      (df["dmr_state_sciMETv2_3842F_excitatory"] == 0))
]

print(df_filtered.shape)  # 119363 rows, 9 columns
df_filtered.to_csv("/secret/path/amethyst/manuscript_analysis_v3/benchmark/allcools/dmr/allcools_sciMETv2_3842F_exc_vs_inh_dmrs.csv", sep="\t", index=False)






