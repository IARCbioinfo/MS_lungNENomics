
"""This module contains functions to compute spatial cross-correlation for spatial transcriptomics samples."""

import pandas as pd
import scanpy as sc
import liana as li
import anndata as ad
import squidpy as sq
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.sparse import csr_matrix, isspmatrix_csr
from tqdm import tqdm


def global_r(x_mat, y_mat, weight):
    #print("Dimensions 2: ", len(x_mat), len(y_mat), weight.shape)
    return ((weight @ x_mat) * y_mat).sum(axis=0)

def zscore(mat, axis=0):
    mat = mat - mat.mean(axis=axis)
    mat = mat / np.sqrt(np.sum(np.power(mat, 2), axis=axis))
    mat = np.clip(mat, -10, 10)

    return np.array(mat)

def spatialdm_weight_norm(weight):
    norm_factor = weight.shape[0] / weight.sum()
    weight = norm_factor * weight
    return weight

def compute_moran(x, y, connectivity):
    x = zscore(x, axis=0)
    y = zscore(y, axis=0)
    weight = spatialdm_weight_norm(connectivity)

    global_stat = global_r(x_mat=x, y_mat=y, weight=weight)

    global_pvals = \
        permutation_pvals(x_mat=x,
                            y_mat=y,
                            weight=weight,
                            global_stat=global_stat,
                            n_perms=10000,
                            seed=1337,
                            )

    return global_stat, global_pvals

def permutation_pvals(x_mat,
                      y_mat,
                      weight,
                      global_stat,
                      n_perms,
                      seed,
                      ):
    rng = np.random.default_rng(seed)

    # initialize mat /w n_perms * number of X->Y
    idx = x_mat.shape[0]
    
    perm_mat = np.zeros((n_perms))

    for perm in tqdm(range(n_perms), disable=True):
        _idx = rng.permutation(idx)
        perm_mat[perm] = global_r(x_mat=x_mat[_idx],
                                        y_mat=y_mat[_idx],
                                        weight=weight)

    global_pvals = ((np.abs(perm_mat) > np.abs(global_stat)).sum(axis=0) + 1) / (n_perms + 1)

    return global_pvals

def compute_spatial_corr(adata_mix, cell_types):

    count_data = adata_mix.to_df(layer="counts_norm")

    connectivity = adata_mix.obsp["spatial_connectivities"]
    connectivity = csr_matrix(connectivity, dtype=np.float32)

    N = len(cell_types)
    
    corr_bv = np.empty((N,N))
    p_bv = np.empty((N,N))

    for i, type1 in enumerate(cell_types):
        print("Step ", i, "/", N)
        for j, type2 in enumerate(cell_types):   
            x = count_data[type1]
            y = count_data[type2]
            #print("Dimensions: ", len(x), len(y), connectivity.shape)
            M_r, p_val = compute_moran(x, y, connectivity)
            corr_bv[i,j] = M_r
            p_bv[i,j] = p_val
                
    return corr_bv, p_bv

def compute_spatial_corr_mixed(adata_mix, cell_types, signal_genes):

    count_data = adata_mix.to_df(layer="counts_norm")

    connectivity = adata_mix.obsp["spatial_connectivities"]
    connectivity = csr_matrix(connectivity, dtype=np.float32)

    N = len(signal_genes)
    M = len(cell_types)
    
    corr_bv = np.empty((N,M))
    p_bv = np.empty((N,M))

    for i, prop in enumerate(signal_genes):
        print("Step ", i, "/", N)
        for j, gene in enumerate(cell_types):   
            x = count_data[prop]
            y = count_data[gene]
            #print("Dimensions: ", len(x), len(y), connectivity.shape)
            M_r, p_val = compute_moran(x, y, connectivity)
            corr_bv[i,j] = M_r
            p_bv[i,j] = p_val
                
    return corr_bv, p_bv

def make_adata_mixed_bv(lung_ad, IRIS_prop, signalling_genes, cell_types):  
    # construct AnnData object with variables of interest
    gene_expr =  lung_ad.to_df(layer="log1p_norm")[signalling_genes]
    cell_prop = IRIS_prop[cell_types]

    expr = pd.concat([gene_expr, cell_prop], axis=1)

    adata_mix = ad.AnnData(csr_matrix(expr, dtype=np.float32), obs = lung_ad.obs)
    adata_mix.obsm = lung_ad.obsm
    adata_mix.uns["spatial"] = lung_ad.uns["spatial"]
    var_names = signalling_genes + list(cell_types)
    adata_mix.var_names = var_names

    adata_mix.layers["counts_norm"] = adata_mix.X

    # Compute neighbor graph
    sq.gr.spatial_neighbors(adata_mix)

    return adata_mix

def compute_spatial_corr_mix(adata_mix, cell_types, signal_genes):

    count_data = adata_mix.to_df(layer="counts_norm")

    connectivity = adata_mix.obsp["spatial_connectivities"]
    connectivity = csr_matrix(connectivity, dtype=np.float32)

    N = len(signal_genes)
    M = len(cell_types)
    
    corr_bv = np.empty((N,M))
    p_bv = np.empty((N,M))

    for i, prop in enumerate(signal_genes):
        print("Step ", i, "/", N)
        for j, gene in enumerate(cell_types):   
            x = count_data[prop]
            y = count_data[gene]
            #print("Dimensions: ", len(x), len(y), connectivity.shape)
            M_r, p_val = compute_moran(x, y, connectivity)
            corr_bv[i,j] = M_r
            p_bv[i,j] = p_val
                
    return corr_bv, p_bv

def plot_correlations(corr_mtx, p_mtx, cell_type_names, gene_names):

    df_corr = pd.DataFrame(corr_mtx, columns=cell_type_names, index=cell_type_names)
    df_pval = pd.DataFrame(p_mtx, columns=cell_type_names, index=cell_type_names)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    sns.heatmap(df_corr, annot = True, cmap="coolwarm", ax=axes[0])
    axes[0].set_title(f"Correlation matrix")
    axes[0].tick_params(axis='y', labelrotation=0)
    #axes[0].tick_params(axis='x', labelrotation=60)


    sns.heatmap(df_pval, annot = True, cmap="coolwarm", ax=axes[1])
    axes[1].set_title(f"p-values")
    axes[1].tick_params(axis='y', labelrotation=0)

    plt.tight_layout()

def plot_correlations_mix(corr_mtx, p_mtx, cell_type_names, gene_names):

    df_corr = pd.DataFrame(corr_mtx, columns=cell_type_names, index=gene_names)
    df_pval = pd.DataFrame(p_mtx, columns=cell_type_names, index=gene_names)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    sns.heatmap(df_corr, annot = True, cmap="coolwarm", ax=axes[0])
    axes[0].set_title(f"Correlation matrix")
    axes[0].tick_params(axis='y', labelrotation=0)
    #axes[0].tick_params(axis='x', labelrotation=60)


    sns.heatmap(df_pval, annot = True, cmap="coolwarm", ax=axes[1])
    axes[1].set_title(f"p-values")
    axes[1].tick_params(axis='y', labelrotation=0)

    plt.tight_layout()