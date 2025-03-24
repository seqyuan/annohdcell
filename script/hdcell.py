#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import os, sys
import Anndata as ad
from annohdcell import read_data, expand_nuclei, expand_nuclei, bin_to_cell
import click


@click.group()
def main() -> None:
    pass

# ------------------------------------------------------------------------------------
@main.command(name="bin2cell")
@click.option('--square_002um', '-d', required=True,
              help="10X HD 2um square bin path")
@click.option('--spatial', '-s', 
              help="spatial path, binned_outputs/square_002um/spatial")
@click.option('--tif', '-t', required=True,
              help="source image file path")
@click.option('--outdir', '-o', required=True,
              help="output dir, one sample per dir")

def hd_cell_sgement(square_002um: str, spatial: str, source_image_path: str, source_image_path: str, outdir: str) -> None:
    adata = read_data(square_002um, spatial, tif, outdir)
    adata = nuclei_detection()
    adata = expand_nuclei(adata)
    cdata = bin_to_cell(adata)

    cdata.var_names_make_unique()
    cdata = cdata[cdata.obs['bin_count']>0] # min 6 bins
    cdata.uns = {}
    del cdata.obsm['spatial_cropped_150_buffer']
    sc.pp.calculate_qc_metrics(cdata,inplace=True)
    cdata.write_h5ad(f'{outdir}/b2c4rds.h5ad')


# ------------------------------------------------------------------------------------
@main.command(name="h5trans")
@click.option('--b2ch5', '-c', default=None, help="bin2cellh5ad")
@click.option('--h54rds', '-o', default=None, help="h5ad path for trans to rds")

def transh5(b2ch5: str) -> None:
    cdata = sc.read_h5ad(b2ch5)
    cdata.var_names_make_unique()
    cdata = cdata[cdata.obs['bin_count']>0] # min 6 bins

    cdata.uns = {}
    del cdata.obsm['spatial_cropped_150_buffer']
    sc.pp.calculate_qc_metrics(cdata,inplace=True)
    cdata.write_h5ad(h54rds)

# ------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

