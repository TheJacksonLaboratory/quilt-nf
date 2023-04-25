#!/bin/bash
#SBATCH -J standalone_QUILT_SJW
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=1

# chromosome
chr=$1
echo ${chr}

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

# create pos file for QUILT
echo "Creating position file"
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools view ${DO_vcf} --regions ${chr} -m2 -M2 -v snps | \
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > QUILT_${chr}_pos.txt
head QUILT_${chr}_pos.txt

# create haplegendsample files for QUILT
echo "Creating haplegendsample files"
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools view ${DO_vcf} --regions ${chr} -m2 -M2 -v snps | \
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools convert --haplegendsample merged_${chr}

# run QUILT
echo "Running QUILT"
singularity exec ${containerDir}/sjwidmay-stitch_nf-stitch.img Rscript --vanilla ${pipeDir}/bin/stitch/run_quilt_DO.R \
        bamlist.txt \
        QUILT_${chr}_pos.txt \
        ${nFounders} \
        ${chr} \
        merged_${chr}.hap.gz \
        merged_${chr}.samples \
        merged_${chr}.legend.gz

echo "Making allele probabilities from QUILT vcf"
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img tabix -p vcf quilt.${chr}.vcf.gz
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%HD]\n' quilt.${chr}.vcf.gz > quilt.alleleprobs.${chr}.txt
