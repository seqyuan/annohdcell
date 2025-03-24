import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import numpy as np
import os
import sys
import bin2cell as b2c

def mask1_h5ad(adata: ad.AnnData):
    mask = ((adata.obs['array_row'] >= 900) & 
        (adata.obs['array_row'] <= 1100) & 
        (adata.obs['array_col'] >= 1100) & 
        (adata.obs['array_col'] <= 1300))
    return mask

def mask2_h5ad(adata: ad.AnnData):
    mask = ((adata.obs['array_row'] >= 400) & 
        (adata.obs['array_row'] <= 600) & 
        (adata.obs['array_col'] >= 600) & 
        (adata.obs['array_col'] <= 800))
    return mask

def read_data(square_002um: str, spatial: str, tif: str, outdir: str, mpp: float = 0.3):
    outdir = os.path.abspath(outdir)
    os.makedirs(os.path.join(outdir, "stardist"), exist_ok=True)
    
    adata = b2c.read_visium(square_002um, 
                        source_image_path=tif, 
                        spaceranger_image_path=spatial)
    adata.var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)
    
    b2c.scaled_he_image(adata, mpp=mpp, save_path=os.path.join(outdir, "stardist/he.tiff"))
    b2c.destripe(adata)
    return adata

def nuclei_detection(adata: ad.AnnData, mpp: float = 0.3):
    b2c.stardist(image_path="stardist/he.tiff", 
             labels_npz_path="stardist/he.npz", 
             stardist_model="2D_versatile_he", 
             prob_thresh=0.01)
    
    b2c.insert_labels(adata, 
                  labels_npz_path="stardist/he.npz", 
                  basis="spatial", 
                  spatial_key="spatial_cropped_150_buffer",
                  mpp=mpp, 
                  labels_key="labels_he")
    return adata

def expand_nuclei(adata: ad.AnnData, mpp: float = 0.3):
    b2c.expand_labels(adata, 
                  labels_key='labels_he', 
                  expanded_labels_key="labels_he_expanded")
    
    b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5, save_path="stardist/gex.tiff")
    b2c.stardist(image_path="stardist/gex.tiff", 
             labels_npz_path="stardist/gex.npz", 
             stardist_model="2D_versatile_fluo", 
             prob_thresh=0.05, 
             nms_thresh=0.5)
    
    b2c.insert_labels(adata, 
                  labels_npz_path="stardist/gex.npz", 
                  basis="array", 
                  mpp=mpp, 
                  labels_key="labels_gex")
    
    b2c.salvage_secondary_labels(adata, 
                             primary_label="labels_he_expanded", 
                             secondary_label="labels_gex", 
                             labels_key="labels_joint")
    return adata

def bin_to_cell(adata: ad.AnnData):
    cdata = b2c.bin_to_cell(adata, 
                          labels_key="labels_joint", 
                          spatial_keys=["spatial", "spatial_cropped_150_buffer"])
    return cdata
