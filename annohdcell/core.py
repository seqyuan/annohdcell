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

    
    return adata[mask]

def mask2_h5ad(adata: ad.AnnData):
    mask = ((adata.obs['array_row'] >= 5400) & 
        (adata.obs['array_row'] <= 5600) & 
        (adata.obs['array_col'] >= 16000) & 
        (adata.obs['array_col'] <= 18000))
    return adata[mask]

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

def vis_roi(adata: ad.AnnData, outdir: str):
    import matplotlib.patches as patches
    fig, ax, = plt.subplots(figsize=(7, 6))
    sc.pl.spatial(adata, color=[None], size=4,
                ax=ax, show=False,
                img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer")

    rect1 = patches.Rectangle((1100, 900), 200, 200, linewidth=1, edgecolor='r', facecolor='none')
    rect2 = patches.Rectangle((16000, 5400), 200, 200, linewidth=1, edgecolor='r', facecolor='none')

    ax.add_patch(rect1)
    ax.add_patch(rect2)
    fig.savefig(f"{outdir}/ROI.he.png")
    fig.savefig(f"{outdir}/ROI.he.pdf")

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

def bin_to_cell(adata: ad.AnnData, outdir: str):
    cdata = b2c.bin_to_cell(adata, 
                          labels_key="labels_joint", 
                          spatial_keys=["spatial", "spatial_cropped_150_buffer"])
    adata.write_h5ad(f'{outdir}/b2c_2um.h5ad')
    cdata.write_h5ad(f'{outdir}/b2c_cell.h5ad')
    return cdata

def vis_stardist(bdata1: ad.AnnData, bdata2: ad.AnnData, 
        outdir: str, 
        oprefix: str,
        color_by=[None, "n_counts", "n_counts_adjusted"],
        img_key: str="0.3_mpp_150_buffer", 
        basis: str="spatial_cropped_150_buffer"):

    fig, axs = plt.subplots(1, len(color_by), figsize=(5*len(color_by)+2, 4))
    for i, v in enumerate(color_by):
         sc.pl.spatial(bdata1, color=v, ax=axs[i], color_map="OrRd", img_key=img_key, basis=basis, show=False)
    plt.tight_layout()
    fig.savefig(f'{outdir}/RO1.{oprefix}.png')
    fig.savefig(f'{outdir}/RO1.{oprefix}.pdf')

    try:
        fig, axs = plt.subplots(1, len(color_by), figsize=(5*len(color_by)-1, 5))
        for i, v in enumerate(color_by):
            sc.pl.spatial(bdata2, color=v, ax=axs[i], color_map="OrRd", img_key=img_key, basis=basis, show=False)
        plt.tight_layout()
        fig.savefig(f'{outdir}/RO2.{oprefix}.png')
        fig.savefig(f'{outdir}/RO2.{oprefix}.pdf')
    except:
        pass

def vis_nuclei_cell(bdata1: ad.AnnData, bdata2: ad.AnnData, 
        outdir: str, 
        oprefix: str="labels_he",
        color_by=[None, "labels_he"],
        img_key: str="0.3_mpp_150_buffer", 
        basis: str="spatial_cropped_150_buffer"):
    
    bdata1 = bdata1[bdata1.obs[color_by[1]]>0]
    bdata1.obs[color_by[1]] = bdata1.obs[color_by[1]].astype(str)
    bdata2 = bdata2[bdata2.obs[color_by[1]]>0]
    bdata2.obs[color_by[1]] = bdata2.obs[color_by[1]].astype(str)

    for ikey in color_by:
        if ikey != None and ikey in bdata1.obs:
            bdata1.obs[ikey] = bdata1.obs[ikey].astype('category')
            bdata2.obs[ikey] = bdata2.obs[ikey].astype('category')

    c_palets1 = custom_palets(bdata1, key=color_by[1])
    c_palets2 = custom_palets(bdata2, key=color_by[1])
    
    fig, axs = plt.subplots(1, len(color_by), figsize=(4*len(color_by)+2, 4))
    for i, v in enumerate(color_by):
        sc.pl.spatial(bdata1, color=v, ax=axs[i], img_key=img_key, basis=basis, show=False, legend_loc='none', palette=c_palets1)  
    plt.tight_layout()
    fig.savefig(f'{outdir}/RO1.{oprefix}.png')
    fig.savefig(f'{outdir}/RO1.{oprefix}.pdf')
    
    try:
        fig, axs = plt.subplots(1, len(color_by), figsize=(4*len(color_by)+2, 4))
        for i, v in enumerate(color_by):
            sc.pl.spatial(bdata2, color=v, ax=axs[i], img_key=img_key, basis=basis, show=False, legend_loc='none', palette=c_palets2)  
        plt.tight_layout()
        fig.savefig(f'{outdir}/RO1.{oprefix}.png')
        fig.savefig(f'{outdir}/RO1.{oprefix}.pdf')
    except:
        pass

def vis_nuclei_cells_heatmap(bdata1: ad.AnnData, bdata2: ad.AnnData, 
        outdir: str, mpp: float = 0.3, 
        im="he", oprefix="stardist_identify_nuclear",
        basis: str="array", spatial_key: str="spatial_cropped_150_buffer"):    

    fig, axs = plt.subplots(1, 2, figsize=(11, 5) )
    crop = b2c.get_crop(bdata1, basis, spatial_key=spatial_key, mpp=mpp)
    rendered = b2c.view_labels(image_path=f"stardist/{im}.tiff", 
                            labels_npz_path=f"stardist/{im}.npz", 
                            crop=crop
                            )
    axs[0].imshow(rendered)

    try:
        crop = b2c.get_crop(bdata2, basis=basis, spatial_key=spatial_key, mpp=mpp)
        rendered = b2c.view_labels(image_path=f"stardist/{im}.tiff", 
                                labels_npz_path=f"stardist/{im}.npz", 
                                crop=crop
                                )
        axs[1].imshow(rendered)
    except:
        pass
    plt.tight_layout()
    fig.savefig(f'{outdir}/ROI.{oprefix}.pdf')
    fig.savefig(f'{outdir}/ROI.{oprefix}.png')

def cell_visualizations(cdata: ad.AnnData, outdir: str):
    mask1 = mask1_h5ad(cdata) 
    mask2 = mask2_h5ad(cdata) 

    ddata1 = cdata[mask1]
    ddata2 = cdata[mask2]

    sc.pl.spatial(ddata1, color=["bin_count", "total_counts"], 
        img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer", 
        save=f'{outdir}/RO1.cell_bin_count.png')

    try:
        sc.pl.spatial(ddata2, color=["bin_count", "total_counts"], 
            img_key="0.3_mpp_150_buffer", basis="spatial_cropped_150_buffer", 
            save=f'{outdir}/RO2.cell_bin_count.png')
    except:
        pass

def custom_palets(bdata, key: str = 'labels_he'):
    from matplotlib.colors import ListedColormap
    import numpy as np
    # 计算类别数量
    n_categories = len(bdata.obs[key].cat.categories)
    # 生成循环扩展的tab20颜色列表（每20种颜色循环一次）
    base_palette = plt.cm.tab20.colors  # 获取tab20的20种基础颜色
    repeat_times = int(np.ceil(n_categories / 20))  # 计算需要重复的次数
    extended_colors = list(base_palette) * repeat_times  # 扩展颜色列表
    c_palets = extended_colors[:n_categories]
    return c_palets
