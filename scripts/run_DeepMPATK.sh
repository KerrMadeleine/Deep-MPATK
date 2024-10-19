#! /bin/bash


# python packages needed
# numpy, matplotlib, csv, sklearn, os, cmcrameri
# you can either install the packages with conda or make a 
# virtual environment with these packages 

# stagyy binary output into csv files 
echo 'READING STAGYY TO BINARY'
#python read_stagyy_bin.py 0 0 0 # for regular models
#python read_stagyy_bin.py 0 1 0 # for high-time-res models


# csv to plotting scalar fields
echo 'PLOTTING T, RHO, VISC, FIELDS'
python plot_fields_spherical.py 0 0 0


writeBL_RACdata=1
haveBL_RACdata=0
writeraeff=0
echo 'PIPELINE USING HIGHER RES MODELS'
#python pipeline.py $haveBL_RACdata $writeBL_RACdata $writeraeff

writeBL_RACdata=0
haveBL_RACdata=0
writeraeff=1
echo 'PIPELINE USING NORMAL MODELS, WRITING RAEFF and CONV REG'
#python pipeline.py $haveBL_RACdata $writeBL_RACdata $writeraeff

# plotting velocity versus time to get Vrms
echo 'PLOTTING VRMS DATA'
#python read_rprof.py $haveBL_RACdata $writeBL_RACdata $writeraeff

echo 'CONV ONSET CALC USING '
#python conv_onset_pars.py $haveBL_RACdata $writeBL_RACdata $writeraeff

writeBL_RACdata=0
haveBL_RACdata=1
writeraeff=0
echo 'PIPELINE'
#python pipeline.py $haveBL_RACdata $writeBL_RACdata $writeraeff


echo 'SNIPPETS'
#python snippets.py $haveBL_RACdata $writeBL_RACdata $writeraeff
