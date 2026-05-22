```bash
project_name="Acer"
dirs=(
    "1_base-filt-output"
    "2_checkambient-output"
    "3_concated-output"
    "4_normalize-output"
    "5_feature-selection-output"
    "6_dim-reduction-output"
    "7_cluster-output"
    "figures"
    "h5_from_fastq2matrix"
    "metadata"
)
for dir in "${dirs[@]}"; do
    mkdir -p "${project_name}/${dir}"
    touch "${project_name}/${dir}/.gitkeep"
done
```

```bash
ln -s "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/drone_rep1/drone_rep1_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Acer/h5_from_fastq2matrix/drone_rep1_singlet.h5"

ln -s "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/drone_rep2/drone_rep2_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Acer/h5_from_fastq2matrix/drone_rep2_singlet.h5"

ln -s "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/queen_rep1/queen_rep1_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Acer/h5_from_fastq2matrix/queen_rep1_singlet.h5"

ln -s "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/worker_rep1/worker_rep1_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Acer/h5_from_fastq2matrix/worker_rep1_singlet.h5"

ln -s "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/singlet-out/worker_rep2/worker_rep2_singlet.h5" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Acer/h5_from_fastq2matrix/worker_rep2_singlet.h5"
```
```bash
ln -s "/home/liuzhiyu/Projects/neo_caste/fastq2matrix/Acer/ref" "/home/liuzhiyu/Projects/neo_caste/downstream_analysis/Acer/ref"
```