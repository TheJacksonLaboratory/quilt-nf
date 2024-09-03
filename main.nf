#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import workflow
// 20230103: Only WGS workflow functional from original cs-nf-pipelines repo
if (params.workflow == "wgs"){
  include {WGS} from './workflows/wgs'
}
if (params.workflow == "quilt"){
  include {QUILT} from './workflows/quilt'
}

// Conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "wgs"){
    WGS()
    }
  if (params.workflow == "quilt"){
    QUILT()
    }
}
