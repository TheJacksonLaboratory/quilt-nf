#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow pipeline for preparing lcWGS data for haplotype reconstruction using QUILT

// Import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/make_quilt_reference_data.nf"
include {DOWNLOAD_REFERENCE_DATA} from "${projectDir}/modules/utility_modules/download_reference_data"
include {FILTER_TO_STRAINS} from "${projectDir}/modules/bcftools/filter_to_strains"

//keep this in in case I need to revive some of the processes
//include {EXPAND_BED} from "${projectDir}/modules/utility_modules/expand_bed.nf"
//include {PILEUPS_TO_BAM} from "${projectDir}/modules/bedtools/filter_bams_to_coverage"
//include {INDEX_FILTERED_BAM} from "${projectDir}/modules/samtools/index_covered_bam"
//include {GENO_PROBS} from "${projectDir}/modules/stitch/genoprobs"
//include {TRIMMOMATIC_PE} from "${projectDir}/modules/utility_modules/trimmomatic"
//include {QUALITY_STATISTICS} from "${projectDir}/modules/utility_modules/quality_stats"
//include {BOWTIE2} from "${projectDir}/modules/bowtie2/bowtie2"
//include {AGGREGATE_STATS} from "${projectDir}/modules/utility_modules/aggregate_stats_wgs"
//include {WRITE_QTL2_FILES} from "${projectDir}/modules/gatk/write_qtl2files"
//include {GENO_PROBS} from "${projectDir}/modules/gatk/genoprobs"

//include {MAKE_B6_VARIANTS} from "${projectDir}/modules/quilt/make_B6_sanger_variants"
//include {MAKE_QUILT_REFERENCE_FILES} from "${projectDir}/modules/quilt/make_haplegendsample"
//include {MAKE_QUILT_MAP} from "${projectDir}/modules/quilt/make_quilt_map"



// help if needed
if (params.help){
    help()
    exit 0
}

// log params
param_log()

// chromosome channel
chrs = Channel.of(1..19,"X")

// main workflow
workflow MAKE_QUILT_REF_DATA {

    // Download reference genome and Sanger SNPs.
    DOWNLOAD_REFERENCE_DATA()
    
    // Gather strains to parse VCF
    if(params.cross_type == 'do'){
      strains = Channel.of('A_J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ')
    } else if(params.cross_type == 'bxd'){
      strains = Channel.of('DBA_2J')
    } else if(params.cross_type == 'cc'){
      strains = Channel.of('A_J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ')
    } else if(params.cross_type == 'het3'){
      strains = Channel.of('C3H_HeJ,BALB_cByJ,DBA_2J')
    }else {
      strains = Channel.of(params.custom_strains)
    }

    // channel for reference files
    reference_files = DOWNLOAD_REFERENCE_DATA.out.ref_data

    // combine with strain set and chromosome channels
    strains_ref_data = reference_files.flatten()
                                      .concat(strains)
                                      .collect()
                                      .combine(chrs)

    // filter Sanger VCF to chromosome and strains
    FILTER_TO_STRAINS(strains_ref_data)

}