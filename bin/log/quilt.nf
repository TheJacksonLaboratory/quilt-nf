def param_log(){

log.info """
______________________________________________________

                QUILT PARAMETER LOG                                                  
                                                                                                 
--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________________________________
--workflow                                          ${params.workflow}
--align_only                                        ${params.align_only}
--cross_type                                        ${params.cross_type}
--covar_file                                        ${params.covar_file}
--downsample_to_cov                                 ${params.downsample_to_cov}
--run_name			                                ${params.run_name}
--read_type                                         ${params.read_type}
--sample_folder                                     ${params.sample_folder}
--pattern                                           ${params.pattern}
--extension                                         ${params.extension}
--concat_lanes                                      ${params.concat_lanes}
-w                                                  ${workDir}
-c                                                  ${params.config}
--pubdir                                            ${params.pubdir}
--organize_by                                       ${params.organize_by}
--ref_fa                                            ${params.ref_fa}
--ref_fa_indices                                    ${params.ref_fa_indices}
--min_pct_hq_reads                                  ${params.min_pct_hq_reads}
--mismatch_penalty                                  ${params.mismatch_penalty}

Project Directory: ${projectDir}
______________________________________________________________________________

"""

}
