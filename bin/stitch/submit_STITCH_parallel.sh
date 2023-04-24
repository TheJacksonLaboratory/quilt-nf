#!/bin/bash

for chr in {1..19}
do
    sbatch barebones_STITCH.sh ${chr} 8
done