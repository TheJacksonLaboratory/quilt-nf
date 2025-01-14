# quilt-nf: A Nextflow pipeline for haplotype reconstruction using low-coverage whole-genome sequencing data

These workflows support the service offered by JAX Genome Technologies for haplotype reconstruction with low-pass WGS. Please contact [Sam Widmayer](mailto:samuel.widmayer@jax.org)  or [Dan Gatti](mailto:dan.gatti@jax.org) for more information.

JAX users are required to have access to the Sumner cluster, and to have Nextflow installed in their home directory. Any setup for external users will require additional support, and those wishing to share these workflows are encouraged to contact the maintainers of this repository.

This pipeline is implemented using [Nextflow](https://www.nextflow.io/), a scalable, reproducible, and increasingly common language used in the development and maintenance of bioinformatics workflows. The modular nature of the workflow is enabled by software containers, such as [Docker](https://www.docker.com/) and [Singularity](https://sylabs.io/singularity), with all the software requirements for executing each step. Specific combinations and versions of software are specified in each container making analyses perfectly reproducible over time as long as the source data is unchanged.

## Overview

```mermaid
flowchart TD
    p0((Sample))
    p1((R/qtl2 Covariate File))
    p2((Reference Haplotypes))
    p3((quilt bin shuffle values file))
    p4((downsampled coverage values file))

    p5[FASTP]:::process
    p6[FASTQC]:::process
    p7[READ GROUPS]:::process
    p8[bwa-mem]:::process
    p9[PICARD SortSam]:::process
    p10[PICARD MarkDuplicates]:::process
    p11[PICARD CollectAlignmentSummaryMetrics]:::process
    p12[PICARD CollectWgsMetrics]:::process
    p13[MULTIQC]:::process
    p14[samtools depth + coverage]:::process

    note1{--align_only}
    note2{Run QUILT}

    p15[DOWNSAMPLE BAM]:::process
    p16[CREATE BAMLIST]:::process
    p17[RUN QUILT]:::process
    p18[QUILT TO R/QTL2 FILES]:::process
    p19[GENOPROBS]:::process

    o1((MULTIQC report)):::output
    o2((PICARD CollectWgsMetrics file)):::output
    o3((PICARD CollectAlignmentSummaryMetrics file)):::output
    o4((samtools depth coverage file)):::output

    o5((quilt VCF file)):::output
    o6((quilt variant filtering summary file)):::output
    o7((R/qtl2 sample genotypes)):::output
    o8((R/qtl2 founder genotypes)):::output
    o9((R/qtl2 genetic and physical maps)):::output
    o10((R/qtl2 Cross Object)):::output
    o11((R/qtl2 36-state/Genotype Probabilities)):::output
    o12((R/qtl2 8-state/Allele Probabilities)):::output


    p0 --> p5
    p1 --> p17
    p2 --> p17
    p3 --> p17
    p4 --> p15
    p5 --> p6
    p6 --> p7
    p6 --> p13
    p7 --> p8
    p8 --> p9
    p9 --> p10
    p10 --> p11
    p10 --> p12
    p11 --> p13
    p12 --> p13
    p10 --> p14
    
    p13 -..-> note1
    p14 -..-> note1
    subgraph align_only [  ]
        note1 --> o1
        note1 --> o2
        note1 --> o3
        note1 --> o4
    end

    p13 -..-> note2
    p14 -..-> note2
    note2 --> p15
    p15 --> p16
    p16 --> p17
    p17 --> o5
    p17 --> p18
    p18 --> o6
    p18 --> o7
    p18 --> o8
    p18 --> o9
    p18 --> p19
    p19 --> o10
    p19 --> o11
    p19 --> o12

classDef output fill:#99e4ff,stroke:#000000,stroke-width:5px,color:#000000
classDef process fill:#00A2DC,stroke:#000000,stroke-width:2px,color:#000000

```

**Execution:**

On the JAX HPC, from within the `quilt-nf` directory:

``` bash
sbatch run_scripts/quilt_DO.sh [run name]
```
A prospective user can write their own run script using the following template:

``` bash
#!/bin/bash
#SBATCH --mail-user={USER.EMAIL}
#SBATCH --job-name=QUILT-NF
#SBATCH --mail-type=END,FAIL
#SBATCH -p compute
#SBATCH -q batch
#SBATCH -t 72:00:00
#SBATCH --mem=1G
#SBATCH --ntasks=1

cd $SLURM_SUBMIT_DIR

# LOAD NEXTFLOW
module use --append /projects/omics_share/meta/modules
module load nextflow

# RUN PIPELINE
nextflow main.nf \
--workflow quilt \
-profile sumner2 \
--sample_folder '{PATH TO DIRECTORY CONTAINTING FASTQ FILES}' \
--gen_org mouse \
--pubdir '{PATH TO DESIRED RESULTS DIRECTORY' \
--extension 'fastq.gz' \ # this is the typical file extension, but see run_scripts/quilt_DO_ddRADseq.sh for alternative example
--pattern="*_R{1,2}*" \ # see above comment
--library_type "seqwell" \ # see above comment
--run_name $1 \
-w '{PATH TO DESIRED NEXTFLOW WORK DIRECTORY}' \ # on JAX, use /flashscratch/{USER} or /flashscratch/{OTHER}
--downsample_to_cov '{PATH TO .CSV WITH COVERAGE VALUES TO DOWNSAMPLE TO}' \
--bin_shuffling_file '{PATH TO .CSV WITH QUILT BIN SHUFFLE RADIUS VALUES}' \
--cross_type 'do' \
--ref_file_dir '{PATH TO DIRECTORY WITH REFERENCE HAPLOTYPE FILES}' \
--covar_file '{PATH TO R/QTL2 COVAR FILE}' \
--comment "This script will run haplotype inference on DO lcWGS data" \
-resume
```
