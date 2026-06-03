```bash
project_name="species/Acer"
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
)
for dir in "${dirs[@]}"; do
    mkdir -p "${project_name}/${dir}"
    touch "${project_name}/${dir}/.gitkeep"
done
```

```bash
# ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/drone_rep1/drone_rep1_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Acer/data/0_h5_from_fastq2matrix/drone_rep1_singlet.h5"

# ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/drone_rep2/drone_rep2_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Acer/data/0_h5_from_fastq2matrix/drone_rep2_singlet.h5"

ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/queen_rep1/queen_rep1_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Acer/data/0_h5_from_fastq2matrix/queen_rep1_singlet.h5"

ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/worker_rep1/worker_rep1_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Acer/data/0_h5_from_fastq2matrix/worker_rep1_singlet.h5"

ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/worker_rep2/worker_rep2_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Acer/data/0_h5_from_fastq2matrix/worker_rep2_singlet.h5"
```
```bash
ln -sf "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/ref" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis_scVI/species/Acer/ref"
```