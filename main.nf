#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import workflow
if (params.workflow == "make_quilt_reference_data"){
  include {MAKE_QUILT_REF_DATA} from './workflows/make_quilt_reference_data'
}
if (params.workflow == "quilt"){
  include {QUILT} from './workflows/quilt'
}

// Conditional to kick off appropriate workflow
workflow{
  if (params.workflow == "make_quilt_reference_data"){
    MAKE_QUILT_REF_DATA()
    }
  if (params.workflow == "quilt"){
    QUILT()
    }
}
