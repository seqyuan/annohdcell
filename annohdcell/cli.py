import click
import os
from .core import (read_data, nuclei_detection, expand_nuclei, vis_roi, vis_stardist, mask1_h5ad, mask2_h5ad,
                  bin_to_cell, save_visualizations, save_intermediate_h5ad, cell_visualizations)
import scanpy as sc

@click.group()
def main() -> None:
    """Command line interface for annohdcell package."""
    pass

@main.command(name="bin2cell")
@click.option('--square_002um', '-d', required=True,
              help="10X HD 2um square bin path")
@click.option('--spatial', '-s', 
              help="spatial path, binned_outputs/square_002um/spatial")
@click.option('--tif', '-t', required=True,
              help="source image file path")
@click.option('--outdir', '-o', required=True,
              help="output directory")
def hd_cell_segment(square_002um: str, spatial: str, tif: str, outdir: str) -> None:
    """Process HD spatial transcriptomics data from bins to cells."""
    os.makedirs(f"{outdir}/stardist", exist_ok=True)
    os.chdir(outdir)
    #os.makedirs("stardist", exist_ok=True)
    
    adata = read_data(square_002um, spatial, tif, outdir)
    vis_roi(adata, outdir)
    bdata1 = mask1_h5ad(adata)
    bdata2 = mask2_h5ad(adata)
    vis_stardist(bdata1, bdata2, 
        outdir=f"{outdir}/stardist", 
        oprefix="n_counts",
        color_by=[None, "n_counts"])

    adata = nuclei_detection(adata)
    bdata1 = mask1_h5ad(adata) 
    bdata2 = mask2_h5ad(adata) 
    vis_nuclei_cell(bdata1, bdata2, 
        outdir=f"{outdir}/stardist",
        oprefix="labels_he",
        color_by=[None, "labels_he"])

    vis_nuclei_cells_heatmap(bdata1, bdata2, outdir=f"{outdir}/stardist", mpp=0.3, im="he", oprefix="stardist_identify_nuclear")

    adata = expand_nuclei(adata)
    bdata1 = mask1_h5ad(adata) 
    bdata2 = mask2_h5ad(adata) 
    vis_nuclei_cell(bdata1, bdata2, 
        outdir=f"{outdir}/stardist",
        oprefix="labels_he_expanded",
        color_by=[None, "labels_he_expanded"])

    vis_nuclei_cell(bdata1, bdata2, 
        outdir=f"{outdir}/stardist",
        oprefix="labels_gex",
        color_by=[None, "labels_gex"])

    vis_nuclei_cells_heatmap(bdata1, bdata2, outdir=f"{outdir}/stardist", mpp=0.3, im="gex", oprefix="b2c_identify_cell")

    cdata = bin_to_cell(adata)
    cell_visualizations(cdata, f'{outdir}/stardist')

    cdata.var_names_make_unique()
    cdata = cdata[cdata.obs['bin_count'] > 0]  # min 6 bins
    cdata.uns = {}
    if 'spatial_cropped_150_buffer' in cdata.obsm:
        del cdata.obsm['spatial_cropped_150_buffer']
    sc.pp.calculate_qc_metrics(cdata, inplace=True)
    cdata.write_h5ad(f'{outdir}/b2c4rds.h5ad')

@main.command(name="h5trans")
@click.option('--b2ch5', '-c', required=True, help="bin2cell h5ad file")
@click.option('--h54rds', '-o', required=True, help="output h5ad path")
def transh5(b2ch5: str, h54rds: str) -> None:
    """Transform h5ad file for R compatibility."""
    cdata = sc.read_h5ad(b2ch5)
    cdata.var_names_make_unique()
    cdata = cdata[cdata.obs['bin_count'] > 0]  # min 6 bins

    cdata.uns = {}
    if 'spatial_cropped_150_buffer' in cdata.obsm:
        del cdata.obsm['spatial_cropped_150_buffer']
    sc.pp.calculate_qc_metrics(cdata, inplace=True)
    cdata.write_h5ad(h54rds)

if __name__ == '__main__':
    main()
