```bash
project_name="Sheng_SA_2020_Hsal"
dirs=(
    "data"
    "data/0_h5_from_fastq2matrix"
    "data/1_base-filt-output"
    "data/2_checkambient-output"
    "data/3_concated-output"
    "data/4_normalize-output"
    "data/5_feature-selection-output"
    "data/6_dim-reduction-output"
    "data/7_cluster-output"
    "figures"
    "figures/1-1_before-filt"
    "figures/1-2_after-filt"
    "figures/2.3-1_checkambient"
    "figures/6_dim-reduction"
    "figures/7_cluster"
    "metadata"
)
for dir in "${dirs[@]}"; do
    mkdir -p "${project_name}/${dir}"
    touch "${project_name}/${dir}/.gitkeep"
done
```

```bash
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013383/GSM4013383_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013383.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013384/GSM4013384_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013384.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013385/GSM4013385_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013385.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013386/GSM4013386_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013386.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013387/GSM4013387_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013387.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013388/GSM4013388_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013388.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013389/GSM4013389_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013389.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013390/GSM4013390_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013390.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013391/GSM4013391_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013391.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013392/GSM4013392_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013392.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Sheng_SA_2020_Hsal/singlet-out/GSM4013393/GSM4013393_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Sheng_SA_2020_Hsal/data/0_h5_from_fastq2matrix/GSM4013393.h5"
```