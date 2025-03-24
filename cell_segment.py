import matplotlib.pyplot as plt
import scanpy as sc
import Anndata as ad
import numpy as np
import os
import sys
import bin2cell as b2c


def mask1_h5ad(adata: ad.AnnData):
    mask = ((adata.obs['array_row'] >= 900) & 
        (adata.obs['array_row'] <= 1100) & 
        (adata.obs['array_col'] >= 1100) & 
        (adata.obs['array_col'] <= 1300)
       )
       
    return mask

def mask2_h5ad(adata: ad.AnnData):
    mask = ((adata.obs['array_row'] >= 400) & 
        (adata.obs['array_row'] <= 600) & 
        (adata.obs['array_col'] >= 600) & 
        (adata.obs['array_col'] <= 800)
    )
       
    return mask

def read_data(square_002um: str, spatial: str, source_image_path: str, source_image_path: str, outdir: str, mpp: float = 0.3):
    #create directory for stardist input/output files
    #outdir = os.path.join(os.path.extend(outdir), "stardist")
    outdir = os.path.abspath(outdir)
    os.chdir(outdir)
    os.makedirs("stardist", exist_ok=True)
    adata = b2c.read_visium(path, 
                        source_image_path = source_image_path, 
                        spaceranger_image_path = spaceranger_image_path
                       )
    adata.var_names_make_unique()
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)

    #likely to be closer to 0.3 for your data
    
    b2c.scaled_he_image(adata, mpp=mpp, save_path="stardist/he.tiff")
    b2c.destripe(adata)

    mask1 = mask1_h5ad(adata)
    mask2 = mask2_h5ad(adata)

    roi_show(adata, mask=mask1, prefix="roi1")
    roi_show(adata, mask=mask2, prefix="roi2")

    return adata
    
def roi_show(adata: ad.AnnData, mask: tuple = mask1, prefix: str = "roi1"):
    #define a mask to easily pull out this region of the object in the future
    bdata = adata[mask]
    sc.pl.spatial(bdata, color=[None, "n_counts", "n_counts_adjusted"], color_map="OrRd",
                img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer", save=f"stardist/{prefix}.pdf")

def nuclei_detection(adata: ad.AnnData, mpp: float = 0.3):
    #run stardist on the cropped image
    b2c.stardist(image_path="stardist/he.tiff", 
             labels_npz_path="stardist/he.npz", 
             stardist_model="2D_versatile_he", 
             prob_thresh=0.01
            )

    b2c.insert_labels(adata, 
                  labels_npz_path="stardist/he.npz", 
                  basis="spatial", 
                  spatial_key="spatial_cropped_150_buffer",
                  mpp=mpp, 
                  labels_key="labels_he"
                 )

    mask1 = mask1_h5ad(adata)
    mask2 = mask2_h5ad(adata)

    labels_he_show(adata, mask=mask1, prefix="roi1")
    labels_he_show(adata, mask=mask2, prefix="roi2")
    return adata
    

def labels_he_show(adata: ad.AnnData, mask: tuple = mask1, prefix: str = "roi1"):
    #define a mask to easily pull out this region of the object in the future
    bdata = adata[mask]

    #the labels obs are integers, 0 means unassigned
    bdata = bdata[bdata.obs['labels_he']>0]
    bdata.obs['labels_he'] = bdata.obs['labels_he'].astype(str)
    sc.pl.spatial(bdata, color=[None, "labels_he"], img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer", 
        save=f"stardist/{prefix}.labels_he.pdf", legend_loc='none')
   
    #the label viewer wants a crop of the processed image
    #get the corresponding coordinates spanning the subset object
    crop = b2c.get_crop(bdata, basis="spatial", spatial_key="spatial_cropped_150_buffer", mpp=mpp)

    rendered = b2c.view_labels(image_path="stardist/he.tiff", 
                            labels_npz_path="stardist/he.npz", 
                            crop=crop
                            )
    plt.imshow(rendered)
    plt.savefig(f'stardist/{prefix}.stardist_identify_nuclear.pdf')
    plt.savefig(f'stardist/{prefix}.stardist_identify_nuclear.png')
   
def expand_nuclei(adata: ad.AnnData, mpp: float = 0.3):
    b2c.expand_labels(adata, 
                  labels_key='labels_he', 
                  expanded_labels_key="labels_he_expanded"
                 )
    labels_he_expand_show(adata, mask=mask1, prefix="roi1")
    labels_he_expand_show(adata, mask=mask2, prefix="roi2")

    b2c.grid_image(adata, "n_counts_adjusted", mpp=mpp, sigma=5, save_path="stardist/gex.tiff")
    b2c.stardist(image_path="stardist/gex.tiff", 
             labels_npz_path="stardist/gex.npz", 
             stardist_model="2D_versatile_fluo", 
             prob_thresh=0.05, 
             nms_thresh=0.5
            )

    b2c.insert_labels(adata, 
                  labels_npz_path="stardist/gex.npz", 
                  basis="array", 
                  mpp=mpp, 
                  labels_key="labels_gex"
                 )
    
    label_gex_show(adata, mask=mask1, prefix="roi1")
    label_gex_show(adata, mask=mask2, prefix="roi2")

    b2c.salvage_secondary_labels(adata, 
                             primary_label="labels_he_expanded", 
                             secondary_label="labels_gex", 
                             labels_key="labels_joint"
                            )
    
    mask1 = mask1_h5ad(adata)
    mask2 = mask2_h5ad(adata)
    label_joint_show(adata, mask=mask1, prefix="roi1")
    label_joint_show(adata, mask=mask2, prefix="roi2")
    
    return adata


def labels_he_expand_show(adata: ad.AnnData, mask: tuple = mask1, prefix: str = "roi1"):
    bdata = adata[mask]
    #the labels obs are integers, 0 means unassigned
    bdata = bdata[bdata.obs['labels_he_expanded']>0]
    bdata.obs['labels_he_expanded'] = bdata.obs['labels_he_expanded'].astype(str)

    sc.pl.spatial(bdata, color=[None, "labels_he_expanded"], img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer", 
        save=f"stardist/{prefix}.labels_he_expand.pdf", legend_loc='none')


def label_gex_show(adata: ad.AnnData, mask: tuple = mask1, prefix: str = "roi1"):
    bdata = adata[mask]

    #the labels obs are integers, 0 means unassigned
    bdata = bdata[bdata.obs['labels_gex']>0]
    bdata.obs['labels_gex'] = bdata.obs['labels_gex'].astype(str)

    sc.pl.spatial(bdata, color=[None, "labels_gex"], img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer", 
        save=f"{prefix}.labels_he_expand.pdf", legend_loc='none')
    crop = b2c.get_crop(bdata, basis="array", mpp=mpp)

    #GEX pops better with percentile normalisation performed
    rendered = b2c.view_labels(image_path="stardist/gex.tiff", 
                            labels_npz_path="stardist/gex.npz", 
                            crop=crop,
                            stardist_normalize=True
                            )
    plt.imshow(rendered)
    plt.savefig(f'stardist/{prefix}.gex_identify_cell.pdf')
    plt.savefig(f'stardist/{prefix}.gex_identify_cell.png')


def label_joint_show(aadata: ad.AnnData, mask: tuple = mask1, prefix: str = "roi1"):
    bdata = adata[mask]

    #the labels obs are integers, 0 means unassigned
    bdata = bdata[bdata.obs['labels_joint']>0]
    bdata.obs['labels_joint'] = bdata.obs['labels_joint'].astype(str)

    sc.pl.spatial(bdata, color=[None, "labels_joint_source", "labels_joint"], 
                img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer", 
                save=f"stardist/{prefix}.labels_he_expand.pdf", legend_loc='none')


def bin_to_cell(adata: ad.AnnData):
    cdata = b2c.bin_to_cell(adata, labels_key="labels_joint", spatial_keys=["spatial", "spatial_cropped_150_buffer"])
    cdata.write_h5ad(f"b2c.h5ad")
    adata.write_h5ad(f"2um.h5ad")

    bin_count_stat_show(cdata)
    return cdata


def bin_count_stat_show(cdata: ad.AnnData):
    mask1 = mask1_h5ad(cdata)
    mask2 = mask2_h5ad(cdata)

    ddata = cdata[mask1]
    sc.pl.spatial(ddata, color=["bin_count", "labels_joint_source"], 
              img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer",
              save=f"stardist/roi1.cell.bin_count.pdf")

    sc.pl.spatial(ddata, color=["total_counts"], size=1.3,
              img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer",
              save=f"stardist/roi1.cell.total_counts.pdf")

    ddata = cdata[mask2]
    sc.pl.spatial(ddata, color=["bin_count", "labels_joint_source"], 
              img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer",
              save=f"stardist/roi2.cell.bin_count.pdf")

    sc.pl.spatial(ddata, color=["total_counts"], size=1.3,
              img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer",
              save=f"stardist/roi2.cell.total_counts.pdf")

    import matplotlib.pyplot as plt
    fig, ax, = plt.subplots(figsize=(7, 6))
    sc.pl.spatial(cdata, color=["total_counts"], size=4,
                ax=ax,
                img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer")

    plt.savefig(f"stardist/total_counts_slice.pdf")

    fig, ax, = plt.subplots(figsize=(4, 4))
    ax = sc.pl.violin(
        adata,
        ["bin_count"],
        jitter=0,
        ax=ax,
        show=False
    )
    ax.text(-0.5, 150, f"cell number: {adata.obs.shape[0]}")
    plt.tight_layout()
    fig.savefig(f"stardist/bin_count_vlnplot.pdf")
    fig.savefig(f"stardist/bin_count_vlnplot.png")

