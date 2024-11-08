#!/bin/bash

for i in ../3_mapping/*.bam
do
    echo "Running QualiMap on $i ..."
    bash ../new_scripts/qualimap.sh $i
done
