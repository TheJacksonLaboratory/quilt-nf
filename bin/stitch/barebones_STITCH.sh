#!/bin/bash
#SBATCH -J standalone_stitch_SJW
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=1

# chromosome
chr=$1

# number of founders
nFounders=$2

# pipeline directory
pipeDir=/projects/compsci/vmp/USERS/widmas/stitch-nf

# project directory; where the original sequence lives
projectDir=${pipeDir}/data/DO_ddRADseq_NovaSeq

# where the bams live as they are spit out from the pipeline
bamDir=${projectDir}/bams

# DO reference variants
DO_vcf=${pipeDir}/bin/stitch/DO_founders.vcf.gz

# singularity directory
containerDir=${pipeDir}/singularity_cache

# create bamlist
ls ${bamDir}/*sorted.marked4_dedup.bam > bamlist.txt

# create pos file for STITCH
if [ -f STITCH_${chr}_pos.txt ]; then
    echo "Position file already exists; skipping"
else
    echo "Creating position file"
    singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools view ${DO_vcf} --regions ${chr} -m2 -M2 -v snps | \
    singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > STITCH_${chr}_pos.txt
fi

# create haplegendsample files for STITCH
if [ -f merged_${chr}.hap.gz ]; then
    echo "Hap files already exist; skipping"
else
    echo "Creating haplegendsample files"
    singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools view ${DO_vcf} --regions ${chr} -m2 -M2 -v snps | \
    singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools convert --haplegendsample merged_${chr}
fi
# run STITCH
if [ -f stitch.${chr}.vcf.gz ]; then
    echo "Already ran STITCH; skipping"
else
    echo "Running STITCH"
    singularity exec ${containerDir}/sjwidmay-stitch_nf-stitch.img Rscript --vanilla ${pipeDir}/bin/stitch/run_stitch_DO.R \
        bamlist.txt \
        STITCH_${chr}_pos.txt \
        ${nFounders} \
        ${chr} \
        merged_${chr}.hap.gz \
        merged_${chr}.samples \
        merged_${chr}.legend.gz
fi

echo "Making allele probabilities from STITCH vcf"
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img tabix -p vcf stitch.${chr}.vcf.gz
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%HD]\n' stitch.${chr}.vcf.gz > stitch.alleleprobs.${chr}.txt