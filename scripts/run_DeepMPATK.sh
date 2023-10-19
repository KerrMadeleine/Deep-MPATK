#! /bin/bash


# python packages needed
# numpy, matplotlib, csv, sklearn, os, cmcrameri
# you can either install the packages with conda or make a 
# virtual environment with these packages 

# stagyy binary output into csv files 
echo 'READING STAGYY TO BINARY'
#python3 read_stagyy_bin.py


# csv to plotting scalar fields
echo 'PLOTTING T, RHO, VISC, FIELDS'
#python3 plot_fields_spherical.py 0 0 0


writeBL_RACdata=1
haveBL_RACdata=0
writeraeff=0
echo 'PIPELINE USING HIGHER RES MODELS'
#python3 pipeline.py $haveBL_RACdata $writeBL_RACdata $writeraeff

writeBL_RACdata=0
haveBL_RACdata=0
writeraeff=1
echo 'PIPELINE USING NORMAL MODELS, WRITING RAEFF and CONV REG'
#python3 pipeline.py $haveBL_RACdata $writeBL_RACdata $writeraeff

# plotting velocity versus time to get Vrms
echo 'PLOTTING VRMS DATA'
#python3 read_rprof.py $haveBL_RACdata $writeBL_RACdata $writeraeff

echo 'CONV ONSET CALC USING '
#python3 conv_onset_pars.py $haveBL_RACdata $writeBL_RACdata $writeraeff

writeBL_RACdata=0
haveBL_RACdata=1
writeraeff=0
echo 'PIPELINE'
#python3 pipeline.py $haveBL_RACdata $writeBL_RACdata $writeraeff


echo 'SNIPPETS'
python3 snippets.py $haveBL_RACdata $writeBL_RACdata $writeraeff
