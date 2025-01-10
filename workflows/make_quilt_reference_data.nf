#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Nextflow pipeline for preparing lcWGS data for haplotype reconstruction using QUILT

// Import modules
include {help} from "${projectDir}/bin/help/wgs.nf"
include {param_log} from "${projectDir}/bin/log/make_quilt_reference_data.nf"
include {DOWNLOAD_INDEX_REFERENCE_DATA} from "${projectDir}/modules/utility_modules/download_reference_data"
include {FILTER_TO_STRAINS} from "${projectDir}/modules/bcftools/filter_to_strains"
include {MAKE_B6_GENOS} from "${projectDir}/modules/make_ref_data/make_B6_genos"
include {MERGE_B6_VCF} from "${projectDir}/modules/bcftools/merge_B6_vcf"
include {PHASE_FOUNDER_VCF} from "${projectDir}/modules/bcftools/phase_founder_vcf"
include {MAKE_REF_HAPLOTYPES} from "${projectDir}/modules/bcftools/make_haplegendsample"
include {MAKE_QUILT_MAP} from "${projectDir}/modules/make_ref_data/make_quilt_map"

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
    DOWNLOAD_INDEX_REFERENCE_DATA()
    
    // Gather strains to parse VCFs
    if(params.cross_type == 'do'){
      strains = Channel.of('A_J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ')
      final_strain_order = Channel.of('A_J,C57BL_6J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ')
    } else if(params.cross_type == 'bxd'){
      strains = Channel.of('DBA_2J')
      final_strain_order = Channel.of('C57BL_6J,DBA_2J')
    } else if(params.cross_type == 'cc'){
      strains = Channel.of('A_J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ')
      final_strain_order = Channel.of('A_J,C57BL_6J,129S1_SvImJ,NOD_ShiLtJ,NZO_HlLtJ,CAST_EiJ,PWK_PhJ,WSB_EiJ')
    } else if(params.cross_type == 'het3'){
      strains = Channel.of('BALB_cByJ,C3H_HeJ,DBA_2J')
      final_strain_order = Channel.of('BALB_cByJ,C57BL_6J,C3H_HeJ,DBA_2J')
    }else {
      strains = Channel.of(params.custom_strains)
      final_strain_order = Channel.of(params.custom_final_strain_order)
    }
    strain_info_channel = strains.concat(final_strain_order).collect()

    // channel for reference files
    reference_files = DOWNLOAD_INDEX_REFERENCE_DATA.out.ref_data

    // combine with strain set and chromosome channels
    strains_ref_data = reference_files.flatten()
                                      .concat(strain_info_channel)
                                      .collect()
                                      .combine(chrs)

    // filter Sanger VCF to chromosome and strains
    FILTER_TO_STRAINS(strains_ref_data)

    // make equivalent B6 variant calls
    MAKE_B6_GENOS(FILTER_TO_STRAINS.out.filtered_sanger_snps)

    // merge B6 vcf into filtered vcf of other strains
    MERGE_B6_VCF(MAKE_B6_GENOS.out.b6_calls)

    // make phased vcf
    PHASE_FOUNDER_VCF(MERGE_B6_VCF.out.complete_vcf)

    // make reference haplotypes
    MAKE_REF_HAPLOTYPES(PHASE_FOUNDER_VCF.out.phased_vcf)

    // make quilt genetic map
    MAKE_QUILT_MAP(PHASE_FOUNDER_VCF.out.phased_vcf)
}