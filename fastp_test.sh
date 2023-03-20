
#!/bin/bash

#SBATCH --job-name=fastp_test
#SBATCH --mail-user=samuel.widmayer@jax.org
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 01:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1

singularity exec docker://sjwidmay/fastp_nf:fastp /fastp -i /fastscratch/STITCH_outputDir/work/e7/1ffbfb7263b926b9a0c400f77ddab6/plexWell-M24_GT23-00380_CTGGTCGT-ATCTATTT_S36_L001_R1_001.fastq.gz_filtered_trimmed -I /fastscratch/STITCH_outputDir/work/e7/1ffbfb7263b926b9a0c400f77ddab6/plexWell-M24_GT23-00380_CTGGTCGT-ATCTATTT_S36_L001_R2_001.fastq.gz_filtered_trimmed -o test.R1.fq -O test.R2.fq

# detect adapters automatically
#        --detect_adapter_for_pe \
# detect and remove duplicates?
#       -D \
# polyG trimming
#        -g \
# base correction analysis
#        -c \
# overrepresentation analysis
#        -p   