# %%

import logging
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import seaborn as sns

# %%

F1_adata = sc.read_h5ad("./Zhang_iScience_2022_Amel/h5_QC/ft-Nm_F1.h5ad")
