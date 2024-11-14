#!/bin/bash
#SBATCH -J 4WC_QUILT_TEST
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 01:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1

# chromosome
chr=$1
echo ${chr}

# pipeline directory
pipeDir=/projects/compsci/vmp/USERS/widmas/stitch-nf

# project directory; where the original sequence lives
projectDir=${pipeDir}/data/4WC_seqwell_NovaSeq_full

# bam directory - doesn't use downsampling
bamDir=${projectDir}/bams

# reference file directory; where each population's 
# chromosome-specific reference files live
refDir=/projects/compsci/vmp/lcgbs_ssif/data/4wc_founders

# singularity directory
containerDir=${pipeDir}/singularity_cache

# QUILT container
QUILT=${containerDir}/sjwidmay-stitch_nf-QUILT.img

# R container
R=${containerDir}/sjwidmay-lcgbs_hr-qtl2_et_al.img

# VariantAnnotation container
VariantAnnotation=${containerDir}/sjwidmay-lcgbs_hr-variantannotation.img

# PLINK and bcftools container
PLINK=${containerDir}/sjwidmay-quilt_nf-plink_bcftools.img

# cross covar file
covar_file=/projects/compsci/vmp/lcgbs_ssif/data/GigaMUGA/4WC_covar.csv

# output directory
outDir=/fastscratch/widmas


### RUN ###

# create bamlist
ls ${bamDir}/*sorted.marked4_dedup.bam > bamlist.txt

echo "Running QUILT"
# run QUILT
singularity exec ${QUILT} Rscript --vanilla ${pipeDir}/bin/quilt/run_quilt.R \
        bamlist.txt \
        ${chr} \
        ${refDir}/chr${chr}.hap.gz \
        ${refDir}/chr${chr}.samples \
        ${refDir}/chr${chr}.legend.gz \
	${covar_file}
        
echo "Pruning markers in tight LD"
# find markers in LD
singularity exec ${PLINK} plink --vcf quilt.${chr}.5000000.7000000.vcf.gz \
    --double-id \
    --maf 0.01 \
    --geno \
    --set-missing-var-ids @:# \
    --chr ${chr} \
    --indep-pairwise 50 20 0.9 \
    --make-bed \
    --out ${outDir}/LDpruned_chr${chr}

# attempt to get a tab separated marker list for bcftools to read
awk '{sub(/:/, "\t"); print}' ${outDir}/LDpruned_chr${chr}.prune.in  > ${outDir}/marker.list.txt

# perform LD pruning
singularity exec ${PLINK} bcftools view quilt.${chr}.5000000.7000000.vcf.gz \
    -R ${outDir}/marker.list.txt \
    --output-type z \
    -o ${outDir}/pruned.quilt.${chr}.vcf.gz

# make qtl2
echo "Writing qtl2 files"
singularity exec ${VariantAnnotation} Rscript --vanilla ${pipeDir}/bin/quilt/prepare_4WC_qtl2_files.R \
    ${refDir}/chr${chr}_phased_snps.vcf.gz \
    ${outDir}/pruned.quilt.${chr}.vcf.gz \
    ${covar_file} \
    ${refDir}/chr${chr}_gen_map.txt \
    ${chr}
    
# make qtl2
echo "Calculating genotype probabilities"
singularity exec ${R} Rscript --vanilla ${pipeDir}/bin/quilt/genoprobs_genail4.R \
    ${chr} \
    chr${chr}_sample_geno.csv \
    chr${chr}_founder_geno.csv \
    chr${chr}_pmap.csv \
    chr${chr}_gmap.csv \
    covar.csv





