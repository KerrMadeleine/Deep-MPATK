import numpy as np
import stagyypythonmodule as spm
import set_DeepMPATK_pars as par
import csv

'''
This script reads in the txt file that supplied information about the onset of convection.
It then writes this information to the pipeline parameter file for use in the next 
stage of the analysis which is to detect plumes.
'''

MODS=par.STAGYY_OUTPUT_MODS

BL_X_KM=[]
BL_X_PT=[]
RAC=[]
N_ONSET=[]
TAU=[]


for i in range(len(MODS)):
    mod0=par.STAGYY_OUTPUT_MODS_0
    mod1=par.STAGYY_OUTPUT_MODS_1
    if i>=mod0 and i<mod1:
        Tablemesh=[]
        with open(par.CONV_ONSET_FILE+MODS[i][:-1]+'.txt') as f:
            reader = csv.reader(f)
            for col in reader:
                Tablemesh.append(col)

        bl_km=float(Tablemesh[-4][1])
        bl_pt=float(Tablemesh[-4][2])
        rac=float(Tablemesh[-4][3])
        n_on=float(Tablemesh[-3][0])
        h_conv=float(Tablemesh[-2][1])
        v_conv=float(Tablemesh[-1][1])

        v_km_per_myr= v_conv*31536*10**6
        tau = h_conv/v_km_per_myr
        
        BL_X_KM.append(bl_km)
        BL_X_PT.append(bl_pt)
        RAC.append(rac)
        N_ONSET.append(n_on)
        TAU.append(tau)
        #print(bl_km,bl_pt,rac,n_on,h_conv,v_conv,tau)
    else:
        BL_X_KM.append(0)
        BL_X_PT.append(0)
        RAC.append(0)
        N_ONSET.append(0)
        TAU.append(0)

pipeline1_par_file = open(par.CONV_ONSET_MOD_PARS+'P1_pars.txt', 'w')
pipeline1_par_file.write(str(BL_X_KM)+'\n')
pipeline1_par_file.write(str(BL_X_PT)+'\n')
pipeline1_par_file.write(str(RAC)+'\n')
pipeline1_par_file.write(str(N_ONSET)+'\n')
pipeline1_par_file.write(str(TAU))
pipeline1_par_file.close()