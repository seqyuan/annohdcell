"""annohdcell package for spatial transcriptomics data processing.

This package provides tools for:
- Reading and processing Visium data
- Nuclei detection and segmentation
- Cell segmentation and bin-to-cell conversion
- Visualization functions
"""

from .core import (
    read_data,
    nuclei_detection,
    expand_nuclei,
    bin_to_cell,
    mask1_h5ad,
    mask2_h5ad
)

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"
