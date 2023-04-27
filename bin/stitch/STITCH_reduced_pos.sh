#!/bin/bash
#SBATCH -J standalone_stitch_SJW
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 20:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=1

# chromosome
chr=$1
echo ${chr}

# number of founders
nFounders=$2

# run description
run_name=$3

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
echo "Creating position file"
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools view --genotype ^het --samples A,B,C,D,E,F,G,H --regions ${chr} -m2 -M2 -v snps --min-ac 2 ${DO_vcf} | \
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > red_STITCH_${chr}_pos.txt
head red_STITCH_${chr}_pos.txt

# create haplegendsample files for STITCH
echo "Creating haplegendsample files"
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools view --genotype ^het --samples A,B,C,D,E,F,G,H --regions ${chr} -m2 -M2 -v snps --min-ac 2 ${DO_vcf}  | \
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools convert --haplegendsample red_merged_${chr}

# run STITCH
echo "Running STITCH"
singularity exec ${containerDir}/sjwidmay-stitch_nf-stitch.img Rscript --vanilla ${pipeDir}/bin/stitch/run_stitch_DO.R \
        bamlist.txt \
        red_STITCH_${chr}_pos.txt \
        ${nFounders} \
        ${chr} \
        red_merged_${chr}.hap.gz \
        red_merged_${chr}.samples \
        red_merged_${chr}.legend.gz \
        stitch_chr${chr}_${run_name}.vcf.gz

echo "Making allele probabilities from STITCH vcf"
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img tabix -p vcf stitch_chr${chr}_${run_name}.vcf.gz
singularity exec ${containerDir}/quay.io-biocontainers-bcftools-1.15--h0ea216a_2.img bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%HD]\n' stitch_chr${chr}_${run_name}.vcf.gz > stitch.alleleprobs.${chr}_${run_name}.txt

echo "Making visualizations and calculating similarity to GigaMUGA"
singularity exec ${containerDir}/lcgbs_hr_qtl2_et_al.sif Rscript /projects/compsci/vmp/USERS/widmas/lcGBS_wf/scripts/stitch_allele_probs.R ${chr} stitch.alleleprobs.${chr}_${run_name}.txt ${run_name}
