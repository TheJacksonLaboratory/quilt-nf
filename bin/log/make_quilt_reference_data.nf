def param_log(){

log.info """
______________________________________________________________________________

                    QUILT REFERENCE DATA PARAMETER LOG                                                  
                                                                                                 
--comment: ${params.comment}

Results Published to: ${params.pubdir}
______________________________________________________________________________
--workflow                                          ${params.workflow}
--cross_type                                        ${params.cross_type}
-w                                                  ${workDir}
-c                                                  ${params.config}
--pubdir                                            ${params.pubdir}

Project Directory: ${projectDir}
______________________________________________________________________________

"""

}
