# %%

import scanpy as sc

# %%

GSM4013383 = sc.read_10x_h5("../Sheng_SA_2020_Hsal/h5_from_fastq2matrix/GSM4013383.h5")

# %%

GSM4013383_doublet = sc.read_10x_h5("../../fastq2matrix/Sheng_SA_2020_Hsal/cellbender-out/GSM4013383/GSM4013383_filtered.h5")
GSM4013383_doublet

# %%

F1 = sc.read_10x_h5("../Zhang_iScience_2022_Amel/h5_from_fastq2matrix/F1.h5")
F1
# %%

F1_doublet = sc.read_10x_h5("../../fastq2matrix/Zhang_iScience_2022_Amel/cellbender-out/GSM5590453/GSM5590453_filtered.h5")
F1_doublet

# %%
