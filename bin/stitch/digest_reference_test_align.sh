#!/bin/bash

##############################################################################################
# Digest mouse genome build 39 with PstI and NlaIII using ddRADseqTools and align test reads
# to the digested genome.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20230508
##############################################################################################

##### SET PARAMS #####

# containers
ddRST_container=/projects/compsci/vmp/lcgbs_ssif/singularity/ddRADseqTools.sif
samtools_container=/projects/compsci/vmp/lcgbs_ssif/singularity/samtools_1.10.sif
bwa_container=/projects/compsci/vmp/USERS/widmas/stitch-nf/singularity_cache/quay.io-biocontainers-bwakit-0.7.17.dev1--hdfd78af_1.img
picard_container=/projects/compsci/vmp/USERS/widmas/stitch-nf/singularity_cache/quay.io-biocontainers-picard-2.26.10--hdfd78af_0.img

# reference genome to be digested
GENFILE=/projects/omics_share/mouse/GRCm39/genome/sequence/ensembl/v105/Mus_musculus.GRCm39.dna.toplevel.fa

# test .fastqs to be aligned
FQ1=/fastscratch/STITCH_outputDir/work/d7/49bf3d2da80a37392a3129dee18aee/ddRADseq-1314_GT23-00337_AAGGA-ACTTGA_S89_R1_filtered_trimmed.fq
FQ2=/fastscratch/STITCH_outputDir/work/d7/49bf3d2da80a37392a3129dee18aee/ddRADseq-1314_GT23-00337_AAGGA-ACTTGA_S89_R2_filtered_trimmed.fq

# output directory
OUTDIR=/fastscratch/widmas

# output digested reference genome
OUTFRAGS=${OUTDIR}/frags_file.fasta

# output digested reference fragment statistics
OUTFRAGSTAT=${OUTDIR}/frags_statistics.csv


##### RUN #####

# digest the reference genome
singularity exec ${ddRST_container} rsitesearch.py --genfile=${GENFILE} --fragsfile=${OUTFRAGS} --fragstfile=${OUTFRAGSTAT} --enzyme1=PstI --enzyme2=NlaIII --minfragsize=125 --maxfragsize=275 --rsfile=/ddRADseqTools/Package/restrictionsites.txt

# rename each fasta line so that the reference can be indexed
sed 's/fragment: //g' ${OUTFRAGS} > ${OUTDIR}/renamed_frags_file.fasta

# index the digested reference genome
# singularity exec ${samtools_container} samtools faidx ${OUTDIR}/renamed_frags_file.fasta
singularity exec ${bwa_container} bwa index ${OUTDIR}/renamed_frags_file.fasta -p ${OUTDIR}/renamed_frags_file.fasta -a bwtsw

# align some test reads (note: they have already been run through fastp
rg=$(cat /fastscratch/STITCH_outputDir/work/64/7984da0109d2dbbfc7381942cd7148/ddRADseq-1314_GT23-00337_AAGGA-ACTTGA_S89_read_group.txt)

singularity exec ${bwa_container} bwa mem -R ${rg} -t 8 -B 8 ${OUTDIR}/renamed_frags_file.fasta ${FQ1} ${FQ2} > ${OUTDIR}/test_fragment_aligned.sam

# index and sort sam
singularity exec ${picard_container} picard -Xmx30G SortSam SO=coordinate INPUT=${OUTDIR}/test_fragment_aligned.sam OUTPUT=${OUTDIR}/test_fragment_aligned.bam VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true