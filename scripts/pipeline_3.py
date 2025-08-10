import numpy as np
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import csv
import os
import sys
import stagyypythonmodule as spm
from cmcrameri import cm

import StagModelClass as smc 
from modelparameters import STAGYY_OUTPUT_FOLDER, MESH_DIR, STAGYY_MODS, STAGYY_MODS_YDIM, STAGYY_MODS_ZDIM, STAGYY_MODS_NMAX, STAGYY_MODS_TCMB, STAGYY_MODS_TM
from modelparameters import TSURF, CORE_RADIUS, DOMAIN_THICKNESS, THERM_EXPANS, GRAV_ACCEL,SPECF_HEAT_CONST_PRESS,THERM_DIFFUS,REF_DENSITY,THERM_COND


"""
pipeline:

X- get T along d_onset evolution
X- get V along d_onset evolution

X- Ra_loc *  (t) = rho * alpha * g * (T_CMB - T_d_onset(t)) * d^3/(kappa * V_d_onset(t)) slice evolution

X- Ra_loc* Kellog plots
X- Ra_loc* Kellog plots with Racrit mask
- Ra_loc* t x y_axis array/mesh save to csv file

"""
rprof_n0 = 0 
rprof_n1 = 57

ord_of_mag_diff = 1
eta_max = 1e25
perc_dev_TBL = 0.1


for m in range(len(STAGYY_MODS)):

    str_label = STAGYY_MODS[m][:-1]

    path = STAGYY_OUTPUT_FOLDER + STAGYY_MODS[m]

    this_model = smc.StagYYModel(path, MESH_DIR, STAGYY_MODS_YDIM[m], STAGYY_MODS_ZDIM[m], STAGYY_MODS_NMAX[m], STAGYY_MODS_TCMB[m], STAGYY_MODS_TM[m], TSURF, CORE_RADIUS, DOMAIN_THICKNESS)
    this_model.set_rprof_dictionary(rprof_n0,rprof_n1)
    this_model.setTimeEvoDF(str_label)
    this_model.setPhysPars(THERM_EXPANS, GRAV_ACCEL,SPECF_HEAT_CONST_PRESS,THERM_DIFFUS,REF_DENSITY,THERM_COND)
    this_model_DF = this_model.readTimeEvoDF(str_label)
    this_model.setOnsetPars(str_label)

    #this_model.makeRaloc_star_Kellog()

    #this_model.makeRaloc_star_Kellog_masked(this_model.Ra_crit_hot)

    this_model.saveRaloc_star_Kellog(Ra_crit = this_model.Ra_crit_hot, saveMaskbool=1, showplotbool=1)

    # individual time step plots
    # for ts in range(0,this_model.nmax+1,100):


    #     T_z_onset_hot = this_model.getHorizontalprofile_z('T',ts,this_model.z_onset_hot)
    #     V_z_onset_hot = this_model.getHorizontalprofile_z('V',ts,this_model.z_onset_hot)

    #     Ra_loc_star_hot = this_model.get_Ra_loc_hot_slice(ts, this_model.z_onset_hot)

    #     #this_model.make_T_profile_plot(ts,perc_dev_TBL,str_label)
    #     #this_model.make_V_profile_plot(ts,eta_max,ord_of_mag_diff,str_label)
