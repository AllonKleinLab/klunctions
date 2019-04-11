mkdir -p bam_saturation_logs
module load gcc samtools/1.9
script="/n/groups/klein/sam/saturation/downsample_bam.py"
conda_dir="/n/groups/klein/sam/miniconda_2019Feb"
pyenv="/n/groups/klein/sam/cnmf"

for f in *bam; do
    s=${f%.bam}
    echo $s
    sbatch -p short --job-name $s --mem 6000 --time 1:00:00 -o bam_saturation_logs/${s}.o -e bam_saturation_logs/${s}.o --wrap """ source ${conda_dir}/bin/activate ${pyenv}; samtools view ${f} | python ${script} ${s}.summary.tsv ${s}.read_counts.tsv ${s}.umi_counts.tsv """
done
