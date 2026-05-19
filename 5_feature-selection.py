# %%

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.sparse import issparse
from scipy.io import mmwrite

# %%


# %%


# ---------------------------------------- Harpegnathos venator ----------------------------------------

# %%

project_name = "Sheng_SA_2020_Hsal"
combined_h5ad = sc.read_h5ad(f"./{project_name}/4_normalize-output/concat.h5ad")

# %%


