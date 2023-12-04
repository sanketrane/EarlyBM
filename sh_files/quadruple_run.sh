#!/bin/bash

while getopts m: flag
do
    case "${flag}" in
        m) modelname=${OPTARG};;
    esac
done

for i in 1 2 3 4
    do
      ./stan_models/${modelname} sample num_warmup=300 num_samples=500 data file=datafiles/Brdu_stanfit.Rdump \
      output file=save_csv/${modelname}_${i}.csv &
    done
