# %%
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import median_abs_deviation

import scanpy as sc

# %%

Project_name = "Sheng_SA_2020_Hsal"
path = ""

# %%

# 注意, 接下来是应该逐个质控最后一起合并, 还是应该先合并然后质控?

