```bash
project_name="species/Zhang_iScience_2022_Amel"
dirs=(
    "data"
    "data/0_h5_from_fastq2matrix"
    "data/1_base-filt-output"
    "data/2_checkambient-output"
    "data/3_concated-output"
    "data/4_feature-selection-output"
    "data/5_scVI-output"
    "data/6_cluster-output"
    "figures"
    "figures/1-1_before-filt"
    "figures/1-2_after-filt"
    "figures/2.3-1_checkambient"
    "figures/5_scVI"
    "figures/6_cluster"
    "metadata"
    "model"
)
for dir in "${dirs[@]}"; do
    mkdir -p "${project_name}/${dir}"
    touch "${project_name}/${dir}/.gitkeep"
done
```

```bash
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Zhang_iScience_2022_Amel/singlet-out/GSM5590453/GSM5590453_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Zhang_iScience_2022_Amel/data/0_h5_from_fastq2matrix/F1.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Zhang_iScience_2022_Amel/singlet-out/GSM5590454/GSM5590454_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Zhang_iScience_2022_Amel/data/0_h5_from_fastq2matrix/F2.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Zhang_iScience_2022_Amel/singlet-out/GSM5590455/GSM5590455_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Zhang_iScience_2022_Amel/data/0_h5_from_fastq2matrix/N1.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Zhang_iScience_2022_Amel/singlet-out/GSM5590456/GSM5590456_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Zhang_iScience_2022_Amel/data/0_h5_from_fastq2matrix/N2.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Zhang_iScience_2022_Amel/singlet-out/GSM5590457/GSM5590457_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Zhang_iScience_2022_Amel/data/0_h5_from_fastq2matrix/Q1.h5"
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Zhang_iScience_2022_Amel/singlet-out/GSM5590458/GSM5590458_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Zhang_iScience_2022_Amel/data/0_h5_from_fastq2matrix/Q2.h5"
```

```bash
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Zhang_iScience_2022_Amel/ref" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/Zhang_iScience_2022_Amel/ref"
```