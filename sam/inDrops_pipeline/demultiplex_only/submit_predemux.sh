# yaml file name
yaml="my_analysis.yaml"

# directory for writing job log files
logdir=logs_predemux

# python environment setup
module load gcc/6.2.0
module load bowtie/1.2.2
module load java/jdk-1.8u112
module load rsem/1.3.0

pyenv=/n/groups/klein/sam/pyndrops
myPython=/n/groups/klein/sam/pyndrops/bin/python
source ${pyenv}/bin/activate ${pyenv}

############################
mkdir -p $logdir

sbatch -p medium -c 8 -N 1 --job-name predemux --mem 16000 -t 72:00:00 -o "$logdir"/predemux.out -e "$logdir"/predemux.out --wrap """ ${myPython} predemultiplexv3fastq.py ${yaml} """

