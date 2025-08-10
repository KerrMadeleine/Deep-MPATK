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
from modelparameters import STAGYY_OUTPUT_FOLDER, MESH_DIR, STAGYY_MODS, STAGYY_MODS_YDIM, STAGYY_MODS_ZDIM, STAGYY_MODS_NMAX, STAGYY_MODS_TCMB, STAGYY_MODS_TM,TSURF, CORE_RADIUS, DOMAIN_THICKNESS
from modelparameters import STAGYY_HIRES_OUTPUT_FOLDER, MESH_DIR, STAGYY_MODS_HIRES, STAGYY_MODS_HIRES_YDIM, STAGYY_MODS_HIRES_ZDIM, STAGYY_MODS_HIRES_NMAX, STAGYY_MODS_HIRES_TCMB, STAGYY_MODS_HIRES_TM,TSURF, CORE_RADIUS, DOMAIN_THICKNESS
from modelparameters import THERM_EXPANS,GRAV_ACCEL,THERM_DIFFUS,REF_DENSITY,THERM_COND

rprof_n0 = 0 
rprof_n1 = 57

ord_of_mag_diff = 1
eta_max = 1e25
perc_dev_TBL = 0.1

"""
pipeline:

X- create a new data frame to store evolutions (low res and high res)

X- time array Myrs evolution
X- mean T evolution
X- mean V evolution
X- get SL thick evolution
X- get convecting region depth evolution
X- get T along SL depth evolution
X- get cold and hot TBL thick evolutions
X- get T along cold TBL evolution
X- get T along hot TBL evolution
X- get V along cold TBL evolution
X- get V along hot TBL evolution
X- get onset convection time for each model

X- Ra = rho * alpha * g * ∆T D ^3 / (kappa * visc) evolution
X- Ra_loc = rho * alpha * g * dT_tbl * d^3/(kappa * visc_tbl) evolution


"""


for m in range(len(STAGYY_MODS)):
    # if '1900' not in STAGYY_MODS[m]:
    #     continue

    # uncomment this :
    # ********************************************************************
    #low res, long time evolution
    str_label = STAGYY_MODS[m][:-1]
    path = STAGYY_OUTPUT_FOLDER + STAGYY_MODS[m]
    this_model = smc.StagYYModel(path, MESH_DIR, STAGYY_MODS_YDIM[m], STAGYY_MODS_ZDIM[m], STAGYY_MODS_NMAX[m], STAGYY_MODS_TCMB[m], STAGYY_MODS_TM[m], TSURF, CORE_RADIUS, DOMAIN_THICKNESS)
    print(str_label)
    # ********************************************************************

    # or this:
    # ********************************************************************
    # high res, short time evolution
    # str_label = STAGYY_MODS[m][:-1]
    # path = STAGYY_HIRES_OUTPUT_FOLDER + STAGYY_MODS_HIRES[m]
    # this_model = smc.StagYYModel(path, MESH_DIR, STAGYY_MODS_HIRES_YDIM[m], STAGYY_MODS_HIRES_ZDIM[m], STAGYY_MODS_HIRES_NMAX[m], STAGYY_MODS_HIRES_TCMB[m], STAGYY_MODS_HIRES_TM[m], TSURF, CORE_RADIUS, DOMAIN_THICKNESS)
    # print(str_label)
    # ********************************************************************



    # uncomment this:

    this_model.set_rprof_dictionary(rprof_n0,rprof_n1)
    time_arr = this_model.getTimeMyrs()
    rad_arr = this_model.getZkm()
    #this_model.createNewEmptyTimeEvoDF(str_label) #uncomment once to create the empty csv files
    this_model.setTimeEvoDF(str_label)




    # # time --------------------
    #this_model.addSeriestoTimeEvoDF(str_label,'time (Myrs)',time_arr)

    # # vrms ----------------------
    # wholedomain_vrms_LM = this_model.getVrmsevolution(2,6) # domain divided into eights
    # this_model.addSeriestoTimeEvoDF(str_label,'vrms (mm/yr)',wholedomain_vrms_LM)

    # # Temperature --------------------
    # meanT_mid_mantle_evolution = this_model.getmeanTevolution(2,5)
    # this_model.addSeriestoTimeEvoDF(str_label,'Tmean (K)',meanT_mid_mantle_evolution)

    # # viscosity --------------------------
    # meanV_mid_mantle_evolution = this_model.getmeanVevolution(2,5)
    # this_model.addSeriestoTimeEvoDF(str_label,'Vmean (Pa s)',meanV_mid_mantle_evolution)


    # # SL thickness -------------------------
    # SL_thick_evo, conv_D_evo, T_SL_evo  = this_model.getSLthickevolution(eta_max, ord_of_mag_diff)
    # this_model.addSeriestoTimeEvoDF(str_label,'SL_thick (m)',SL_thick_evo)
    # this_model.addSeriestoTimeEvoDF(str_label,'D_conv (m)',conv_D_evo)
    # this_model.addSeriestoTimeEvoDF(str_label,'T_SL_base (K)',T_SL_evo)


    # # TBL thicknesses ----------------------------
    # hot_TBL_z_evo, cold_TBL_z_evo, hot_TBL_T_evo, hot_TBL_V_evo, cold_TBL_T_evo, cold_TBL_V_evo = this_model.getTBL_z_evolution(perc_dev_TBL, str_label)
    # this_model.addSeriestoTimeEvoDF(str_label,'hot_TBL_z (m)',hot_TBL_z_evo)
    # this_model.addSeriestoTimeEvoDF(str_label,'cold_TBL_z (m)',cold_TBL_z_evo)
    # this_model.addSeriestoTimeEvoDF(str_label,'hot_TBL_T (m)',hot_TBL_T_evo)
    # this_model.addSeriestoTimeEvoDF(str_label,'cold_TBL_T (m)',cold_TBL_T_evo)
    # this_model.addSeriestoTimeEvoDF(str_label,'hot_TBL_V (m)',hot_TBL_V_evo)
    # this_model.addSeriestoTimeEvoDF(str_label,'cold_TBL_V (m)',cold_TBL_V_evo)
    # this_model.writeTimeEvoDF(str_label)


    # # -------------- overturn times -------------------
    # time_overturns = this_model.getTimeOverturns(str_label)
    # this_model.addSeriestoTimeEvoDF(str_label,'time (overturns)',time_overturns)


    # write all this additions to the data frame to the csv file:
    #this_model.writeTimeEvoDF(str_label)



    this_model_DF = this_model.readTimeEvoDF(str_label)

    # # # ------------- Whole mantle Ra -------------------
    # # #rho * alpha * g * ∆T D ^3 / (kappa * visc)
    # Ra_evo = (REF_DENSITY * THERM_EXPANS * GRAV_ACCEL * (this_model.Tcmb - this_model_DF['T_SL_base (K)']) * this_model_DF['D_conv (m)']**3)/(THERM_DIFFUS * this_model_DF['Vmean (Pa s)'])
    # this_model.addSeriestoTimeEvoDF(str_label,'Ra_WM',Ra_evo)

    # # # ------------- Ra_eff -------------------
    # # #rho * alpha * g * ∆T D ^3 / (kappa * visc)
    # Ra_evo = (REF_DENSITY * THERM_EXPANS * GRAV_ACCEL * (this_model.Tcmb - this_model_DF['T_SL_base (K)']) * this_model_DF['D_conv (m)']**3)/(THERM_DIFFUS * this_model_DF['Vmean (Pa s)'])
    # this_model.addSeriestoTimeEvoDF(str_label,'Ra_eff',Ra_evo)

    # # # -------------Local Ra (hot) ---------------------

    # Ra_loc_hot_evo = (REF_DENSITY * THERM_EXPANS * GRAV_ACCEL * (this_model.Tcmb - this_model_DF['hot_TBL_T (m)']) * this_model_DF['hot_TBL_z (m)']**3)/(THERM_DIFFUS * this_model_DF['hot_TBL_V (m)'])
    # this_model.addSeriestoTimeEvoDF(str_label,'Ra_loc_hot',Ra_loc_hot_evo)

    # # # --------------Local Ra (cold) --------------------

    #cold_TBL_thick_evo = (this_model.D - this_model_DF['SL_thick (m)']) - this_model_DF['cold_TBL_z (m)']
    #Ra_loc_cold_evo = (REF_DENSITY * THERM_EXPANS * GRAV_ACCEL * (this_model_DF['cold_TBL_T (m)'] - this_model_DF['T_SL_base (K)']) * cold_TBL_thick_evo **3)/(THERM_DIFFUS * this_model_DF['cold_TBL_V (m)'])
    #this_model.addSeriestoTimeEvoDF(str_label,'Ra_loc_cold',Ra_loc_cold_evo)


    this_model_DF = this_model.readTimeEvoDF(str_label)

    # fig, ax = plt.subplots(1)
    # fig.set_size_inches(9.5,7)
    # ax.plot( time_arr,  this_model_DF['Ra_WM'], color = 'k')
    # ax.plot( time_arr,  this_model_DF['Ra_loc_hot'], color = 'r')
    # ax.plot( time_arr,  this_model_DF['Ra_loc_cold'], color = 'b')
    # ax.set_yscale('log')
    # #ax.set_xlim([0,5])
    # plt.show()
    # plt.clf()
    # plt.close()

    ## --------------------------------------------------------------------------------------
    ## --------------------------------------------------------------------------------------



# plot all models on same axes

# fig, ax = plt.subplots(1)
# fig.set_size_inches(9.5,7)
# for m in range(len(STAGYY_MODS)):
#     str_label = STAGYY_MODS[m][:-1]
#     print(str_label)
#     path = STAGYY_OUTPUT_FOLDER + STAGYY_MODS[m]
#     this_model = smc.StagYYModel(path, MESH_DIR, STAGYY_MODS_YDIM[m], STAGYY_MODS_ZDIM[m], STAGYY_MODS_NMAX[m], STAGYY_MODS_TCMB[m], STAGYY_MODS_TM[m], TSURF, CORE_RADIUS, DOMAIN_THICKNESS)
#     this_model_DF = this_model.readTimeEvoDF(str_label)
#     ax.plot(this_model_DF['time (overturns)'],  this_model_DF['Ra_WM'])
# ax.set_yscale('log')
# #ax.set_xlim([0,5])
# plt.show()
# plt.clf()
# plt.close()






