# %%

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import median_abs_deviation

import scanpy as sc
import scvi

sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=150, facecolor="white", frameon=False)

# %%

F1_adata = sc.read_10x_h5("./Zhang_iScience_2022_Amel/h5_from_fastq2matrix/F1.h5")
F1_adata

# %%

F2_adata = sc.read_10x_h5("./Zhang_iScience_2022_Amel/h5_from_fastq2matrix/F2.h5")
F2_adata

# %%


