#!/bin/bash

################################################################################
# Submits script to digest the reference genome with a list of enzymes supplied 
# in concert with PstI. Note the bin/stitch/digest_reference.sh 
# accepts any two enzymes.
#
# Sam Widmayer
# samuel.widmayer@jax.org
# 20230511
################################################################################

# list of other enzymes
ENZYME_LIST=/projects/compsci/vmp/lcgbs_ssif/data/enzyme_test_digest.txt
DIGEST_SCRIPT=/projects/compsci/vmp/USERS/widmas/stitch-nf/bin/stitch/digest_reference.sh 

while read line; do
  sbatch ${DIGEST_SCRIPT} $line PstI
done < ${ENZYME_LIST}
