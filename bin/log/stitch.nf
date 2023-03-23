def param_log(){
  
log.info """
______________________________________________________

                STITCH PARAMETER LOG

--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________
--workflow                      ${params.workflow}
--gen_org                       ${params.gen_org}
--do_mice                       ${params.do_mice}
--covar_file                    ${params.covar_file}
--seq_method                    ${params.seq_method}
--read_type                     ${params.read_type}
--sample_folder                 ${params.sample_folder}
--pattern                       ${params.pattern}
--extension                     ${params.extension}
--concat_lanes                  ${params.concat_lanes}
-w                              ${workDir}
-c                              ${params.config}
--pubdir                        ${params.pubdir}
--organize_by                   ${params.organize_by}
--ref_fa                        ${params.ref_fa}
--ref_fa_indices                ${params.ref_fa_indices}
--min_pct_hq_reads              ${params.min_pct_hq_reads}
--mismatch_penalty              ${params.mismatch_penalty}

Project Directory: ${projectDir}
______________________________________________________
"""

}
