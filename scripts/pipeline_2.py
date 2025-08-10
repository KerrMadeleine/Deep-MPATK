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

"""
pipeline:

- use the hires models to get: 
- Ra_loc @ onset time = Ra_crit x ntime
- t_onset x ntime
- d_onset x ntime
- put the onset statistics into the long-time, low-res models.

"""

rprof_n0 = 0 
rprof_n1 = 57

ord_of_mag_diff = 1
eta_max = 1e25
perc_dev_TBL = 0.1

for m in range(len(STAGYY_MODS)):
    # if '1300' not in STAGYY_MODS[m]:
    #     continue

    # ********************************************************************
    # high res, short time evolution
    str_label = STAGYY_MODS[m][:-1]
    path = STAGYY_HIRES_OUTPUT_FOLDER + STAGYY_MODS_HIRES[m]
    this_model = smc.StagYYModel(path, MESH_DIR, STAGYY_MODS_HIRES_YDIM[m], STAGYY_MODS_HIRES_ZDIM[m], STAGYY_MODS_HIRES_NMAX[m], STAGYY_MODS_HIRES_TCMB[m], STAGYY_MODS_HIRES_TM[m], TSURF, CORE_RADIUS, DOMAIN_THICKNESS)
    print(str_label)
    # ********************************************************************

    this_model.set_rprof_dictionary(rprof_n0,rprof_n1)
    time_arr = this_model.getTimeMyrs()
    rad_arr = this_model.getZkm()
    #this_model.createNewEmptyTimeEvoDF(str_label) #uncomment once to create the empty csv files
    this_model.setTimeEvoDF(str_label)


    # # time --------------------
    # this_model.addSeriestoTimeEvoDF(str_label,'time (Myrs)',time_arr)

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


    # #------ plot evolutions as created ----------------------------------------------------
    ## ----------------use for high res------------------------------------------------
    this_model_DF = this_model.readTimeEvoDF(str_label)

    # --------- onset time (hot) ********************************************************
    times_Myrs = this_model_DF['time (Myrs)']
    hot_TBL_evo = this_model_DF['hot_TBL_z (m)']
    t_idx_onset_cold = np.argmax(hot_TBL_evo)
    t_onset_hot = times_Myrs[t_idx_onset_cold]
    fig, ax = plt.subplots(1)
    fig.set_size_inches(9.5,7)
    ax.plot(times_Myrs,hot_TBL_evo)
    plt.show()
    plt.clf()
    plt.close()
    # #this_model.addSeriestoTimeEvoDF(str_label,'time_onset_hot (Myrs)',np.array([t_onset_hot]*(this_model.nmax+1)))

    # --------- onset time (cold) ********************************************************
    times_Myrs = this_model_DF['time (Myrs)']
    cold_TBL_evo = this_model_DF['cold_TBL_z (m)']
    t_idx_onset_hot = np.argmin(cold_TBL_evo)
    t_onset_cold = times_Myrs[t_idx_onset_hot]
    # #this_model.addSeriestoTimeEvoDF(str_label,'time_onset_cold (Myrs)',np.array([t_onset_cold]*(this_model.nmax+1)))


    # --------- onset hot TBL thickness ********************************************************
    final_time = times_Myrs[this_model.nmax]
    time_threshold_polyfit = final_time * 0.9
    time_below_threshold = [times_Myrs[i] for i in range(this_model.nmax+1) if times_Myrs[i]<time_threshold_polyfit]
    hotTBL_below_threshold = [this_model_DF['hot_TBL_z (m)'][i] for i in range(this_model.nmax+1) if times_Myrs[i]<time_threshold_polyfit]

    approx_time_coeffs = np.polyfit(hotTBL_below_threshold, time_below_threshold, deg=2)
    approx_time_poly = np.poly1d(approx_time_coeffs)
    approx_time_hot = approx_time_poly(this_model_DF['hot_TBL_z (m)'])

    z_onset_hot = np.interp(t_onset_hot, approx_time_hot, this_model_DF['hot_TBL_z (m)'])
    #this_model.addSeriestoTimeEvoDF(str_label,'hot_TBL_onset (m)', np.array([z_onset_hot]*(this_model.nmax+1)))


    # --------- onset cold TBL thickness ********************************************************
    final_time = times_Myrs[this_model.nmax]
    time_threshold_polyfit = final_time * 0.9
    time_below_threshold = [times_Myrs[i] for i in range(this_model.nmax+1) if times_Myrs[i]<time_threshold_polyfit]
    coldTBL_below_threshold = [this_model_DF['cold_TBL_z (m)'][i] for i in range(this_model.nmax+1) if times_Myrs[i]<time_threshold_polyfit]

    approx_time_coeffs = np.polyfit(coldTBL_below_threshold, time_below_threshold, deg=2)
    approx_time_poly = np.poly1d(approx_time_coeffs)
    approx_time_cold = approx_time_poly(this_model_DF['cold_TBL_z (m)'])

    z_onset_cold = np.interp(t_onset_cold, approx_time_cold, this_model_DF['cold_TBL_z (m)'])
    #this_model.addSeriestoTimeEvoDF(str_label,'cold_TBL_onset (m)', np.array([z_onset_cold]*(this_model.nmax+1)))


    # -------- high res plotting for z_onset -----------********************************************************
    # fig, ax = plt.subplots(1)
    # fig.set_size_inches(9.5,7)
    # ax.plot( times_Myrs,  this_model_DF['cold_TBL_z (m)'], color = 'k')
    # ax.plot( approx_time_cold,  this_model_DF['cold_TBL_z (m)'], color = 'b')
    # ax.plot(np.ones(this_model.nmax+1) * time_threshold_polyfit,this_model_DF['cold_TBL_z (m)'], color='r')
    # ax.plot(np.ones(this_model.nmax+1) * t_onset_cold,this_model_DF['cold_TBL_z (m)'], color='y')
    # ax.plot( times_Myrs,  np.ones(this_model.nmax+1) * z_onset_cold, color = 'g')

    # ax.plot( times_Myrs,  this_model_DF['hot_TBL_z (m)'], color = 'k')
    # ax.plot( approx_time_hot,  this_model_DF['hot_TBL_z (m)'], color = 'b')
    # ax.plot(np.ones(this_model.nmax+1) * time_threshold_polyfit,this_model_DF['hot_TBL_z (m)'], color='r')
    # ax.plot(np.ones(this_model.nmax+1) * t_onset_hot,this_model_DF['hot_TBL_z (m)'], color='y')
    # ax.plot( times_Myrs,  np.ones(this_model.nmax+1) * z_onset_hot, color = 'g')

    # plt.show()
    # plt.clf()
    # plt.close()


    this_model_DF = this_model.readTimeEvoDF(str_label)

    # # -------------Local Ra crit (hot) ---------------------********************************************************
    z_onset_hot = this_model_DF['hot_TBL_onset (m)'][0]
    T_at_z_onset_t_onset_hot = this_model_DF['Tmean (K)'][t_idx_onset_hot] #this_model.getT_zon_ton(z_onset_hot,t_idx_onset_hot)
    V_at_z_onset_t_onset_hot = this_model.getV_zon_ton(z_onset_hot,t_idx_onset_hot)


    Ra_loc_hot_crit = (REF_DENSITY * THERM_EXPANS * GRAV_ACCEL * (this_model.Tcmb - T_at_z_onset_t_onset_hot) * z_onset_hot**3)/(THERM_DIFFUS * V_at_z_onset_t_onset_hot)
    print('Ra crit hot:', Ra_loc_hot_crit)
    this_model.addSeriestoTimeEvoDF(str_label,'Ra_crit_hot',np.array([Ra_loc_hot_crit]*(this_model.nmax+1)))

    # # --------------Local Ra crit (cold) --------------------********************************************************

    z_onset_cold = this_model_DF['cold_TBL_onset (m)'][0]
    T_at_z_onset_t_onset_cold = this_model.getT_zon_ton(z_onset_cold,t_idx_onset_cold)
    V_at_z_onset_t_onset_cold = this_model.getV_zon_ton(z_onset_cold,t_idx_onset_cold)

    SL_thick_tonset = this_model_DF['SL_thick (m)'][t_idx_onset_cold]
    d_coldTBL_tonset = (this_model.D - SL_thick_tonset) - z_onset_cold
    T_SL_base_tonset = this_model_DF['T_SL_base (K)'][t_idx_onset_cold]
    Ra_loc_cold_crit = (REF_DENSITY * THERM_EXPANS * GRAV_ACCEL * (T_at_z_onset_t_onset_cold - T_SL_base_tonset) * d_coldTBL_tonset**3)/(THERM_DIFFUS * V_at_z_onset_t_onset_cold)
    print('Ra crit cold:', Ra_loc_cold_crit)
    this_model.addSeriestoTimeEvoDF(str_label,'Ra_crit_cold',np.array([Ra_loc_cold_crit]*(this_model.nmax+1)))


    this_model.writeTimeEvoDF(str_label)
    this_model_DF = this_model.readTimeEvoDF(str_label)

    #****************************************************************************************************************








    # ----------------- SWAP HIRES to LORES ----------------------------------
    # ---------------------------------------------------------------------------



    str_label = STAGYY_MODS[m][:-1]
    path = STAGYY_OUTPUT_FOLDER + STAGYY_MODS[m]
    this_model = smc.StagYYModel(path, MESH_DIR, STAGYY_MODS_YDIM[m], STAGYY_MODS_ZDIM[m], STAGYY_MODS_NMAX[m], STAGYY_MODS_TCMB[m], STAGYY_MODS_TM[m], TSURF, CORE_RADIUS, DOMAIN_THICKNESS)
    print(str_label)

    this_model.set_rprof_dictionary(rprof_n0,rprof_n1)
    time_arr = this_model.getTimeMyrs()
    rad_arr = this_model.getZkm()
    #this_model.createNewEmptyTimeEvoDF(str_label) #uncomment once to create the empty csv files
    this_model.setTimeEvoDF(str_label)
    this_model_DF = this_model.readTimeEvoDF(str_label)

    times_Myrs = this_model_DF['time (Myrs)']
    hot_TBL_evo = this_model_DF['hot_TBL_z (m)']
    t_idx_onset_cold = np.argmax(hot_TBL_evo)
    t_onset_hot = times_Myrs[t_idx_onset_cold]
    
    # fig, ax = plt.subplots(1)
    # fig.set_size_inches(9.5,7)
    # ax.plot(times_Myrs,hot_TBL_evo)
    # plt.show()
    # plt.clf()
    # plt.close()

    print('this model nmax = ', this_model.nmax)

    # add t_onset and z_onset (hot and cold) to the low res dataframes
    this_model.addSeriestoTimeEvoDF(str_label,'time_onset_hot (Myrs)',np.array([t_onset_hot]*(this_model.nmax+1)))
    this_model.addSeriestoTimeEvoDF(str_label,'time_onset_cold (Myrs)',np.array([t_onset_cold]*(this_model.nmax+1)))
    this_model.addSeriestoTimeEvoDF(str_label,'hot_TBL_onset (m)', np.array([z_onset_hot]*(this_model.nmax+1)))
    this_model.addSeriestoTimeEvoDF(str_label,'cold_TBL_onset (m)', np.array([z_onset_cold]*(this_model.nmax+1)))
    this_model.addSeriestoTimeEvoDF(str_label,'Ra_crit_hot',np.array([Ra_loc_hot_crit]*(this_model.nmax+1)))
    this_model.addSeriestoTimeEvoDF(str_label,'Ra_crit_cold',np.array([Ra_loc_cold_crit]*(this_model.nmax+1)))

    this_model.writeTimeEvoDF(str_label)