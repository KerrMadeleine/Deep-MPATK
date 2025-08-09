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
todo: 
X- use the Kellog_Ralocstar and Racrit to interpolate plume material widths
X- between ts_0 and ts_1, get the width in m of all plume material
X- make clustered plots

X- move clustered kellog plots to module
X- make snippets clustered kellogs and statistics plots (t_0, t_1)
X- make that into a method
X- make a method that extracts k=2 statistics between t0 and t1
X- make a method that averages the Ra between t0 and t1
X- make a method that plots Ra vs k=2 statistics for a number of snippets
X- make a methods that plots Ra vs k=2 stats for all models and numbers of snippet

- between ts_0 and ts_1, get the width in m of all plume material NOT touching the walls

X- power law fit to heads and conduits vs Ra
X- power law width vs Ra plot with Lithgow Bertelloni
- viscosity plots corn syrup vs arrhenius fluid

- update figure 4a
- figure 5
- figure 7
- figure 9

"""
rprof_n0 = 0 
rprof_n1 = 57

ord_of_mag_diff = 1
eta_max = 1e25
perc_dev_TBL = 0.1

k_cluster = 5
k_choice = 2

ts_start = 100
ts_per_snippet = 100


RA_EFF = list()
RA_EFF_HEAD = list()
RA_EFF_COND = list()
HEAD_STATS = list()
CONDUIT_STATS = list()
ALL_HEAD = list()
ALL_CONDUITS = list()


# ------------------- Looping through each model instance ------------------

for m in range(len(STAGYY_MODS)):

    if '1300' in STAGYY_MODS[m]:
        continue 
    if '1400' in STAGYY_MODS[m]:
        continue    

    # string lavbl for accessing the right dataframe 
    str_label = STAGYY_MODS[m][:-1]

    path = STAGYY_OUTPUT_FOLDER + STAGYY_MODS[m]

    this_model = smc.StagYYModel(path, MESH_DIR, STAGYY_MODS_YDIM[m], STAGYY_MODS_ZDIM[m], STAGYY_MODS_NMAX[m], STAGYY_MODS_TCMB[m], STAGYY_MODS_TM[m], TSURF, CORE_RADIUS, DOMAIN_THICKNESS)
    
    # set the rprof dictionary as an object attribute
    this_model.set_rprof_dictionary(rprof_n0,rprof_n1)

    # set the dataframe in the csv file names [str_label].csv as an object attribute
    this_model.setTimeEvoDF(str_label)

    # set other physicsal parameters to be attributes
    this_model.setPhysPars(THERM_EXPANS, GRAV_ACCEL,SPECF_HEAT_CONST_PRESS,THERM_DIFFUS,REF_DENSITY,THERM_COND)
    
    # make a local copy of the data frame
    this_model_DF = this_model.TimeEvoDFs_dict[str_label]

    this_model.setOnsetPars(str_label)
    cmap_rb , cmap_bin = this_model.getColorMaps()

    Ra_loc_star_Kellog = this_model.readKellog_Ralocstar(plotBool=0)
    Ra_loc_star_Mask = this_model.readKellog_Mask(plotBool=0)

    times_Myrs = this_model_DF['time (Myrs)']


    ts_bounds = this_model.create_snippet_bounds(ts_per_snippet, ts_start)


    Ra_eff_snippets = list()
    Ra_eff_head_snippets = list()
    Ra_eff_cond_snippets = list()
    head_stats_snippets = list()
    conduit_stats_snippets = list()
    all_heads = list()
    all_conduits = list()

    # ------------- Looping through each snippet ------------
    for ts_0, ts_1 in ts_bounds:

        print('snippet from: ts= ', ts_0, ' to ts = ', ts_1)

        Ra_eff_t0_t1 = this_model.getRa_bw_t0_t1(str_label,ts_0,ts_1)

        #Between ts_0 and ts_1, find all the plumes total in the snippet and all the plumes per time
        # step in the snippet (int, list, list ,list)
        number_of_all_plumes_snippet, all_plumes_snippet, number_of_each_plumes_slice, each_plume_slice = this_model.getPlumeCounts(ts_0,ts_1, Ra_loc_star_Kellog,showMiniplotsBool=0)

        heads, conduits = this_model.getListOfAllHeads_n_Conds(all_plumes_snippet)

        Ra_eff_head = [Ra_eff_t0_t1]*len(heads)
        Ra_eff_cond = [Ra_eff_t0_t1]*len(conduits)

        # get clusterings and statistics on all plumes in the snippet from k=1 to k=k_cluster
        kmeans_out = this_model.Kmeans_all_plumes( kmax=k_cluster , allplumes=all_plumes_snippet)
        
        # get all centroids from k=1 to k=k_cluster
        Centroids = kmeans_out[0]
        
        plume_conduit_stats, plume_head_stats = this_model.getHeadConduitStats(kmeans_out)
        #this_model.makeStatsPlots(kmeans_out)

        Ra_eff_snippets.append(Ra_eff_t0_t1)
        Ra_eff_head_snippets.append(Ra_eff_head)
        Ra_eff_cond_snippets.append(Ra_eff_cond)
        head_stats_snippets.append(plume_head_stats)
        conduit_stats_snippets.append(plume_conduit_stats)
        all_heads.append(heads)
        all_conduits.append(conduits)

        k_idx = k_cluster-1

        if len(Centroids)<k_cluster:
            print('len(centroids) = ', len(Centroids))
            print('ktest = ', k_cluster)
            print('len(Centroids)<k_cluster')
            continue

        #this_model.plotGroupedKellog(ts_0,ts_1,k_choice,Ra_loc_star_Mask,number_of_each_plumes_slice,each_plume_slice,Centroids)


    # ------------------- FOR ALL SNIPPETS FOR 1 MODEL -------------------

    # single model head and conduit size evolution
    fig, ax = plt.subplots(1)
    fig.set_size_inches(5,3)
    fig.suptitle(str(int(this_model.Tm_0)))
    ax.scatter(Ra_eff_snippets,[snip[0] for snip in head_stats_snippets],color='k',s=[snip[2] for snip in head_stats_snippets], marker='o')
    ax.scatter(Ra_eff_snippets,[snip[0] for snip in conduit_stats_snippets], color='g',s=[snip[2] for snip in conduit_stats_snippets], marker ='o')
    #for sn_idx in range(len(ts_bounds)):
    #    ax.scatter([Ra_eff_snippets[sn_idx]]*len(all_heads[sn_idx]), all_heads[sn_idx], color ='k', alpha = 0.05)
    #    ax.scatter([Ra_eff_snippets[sn_idx]]*len(all_conduits[sn_idx]), all_conduits[sn_idx], color ='g', alpha = 0.05)
    ax.errorbar(Ra_eff_snippets,[snip[0] for snip in head_stats_snippets],yerr=[snip[1] for snip in head_stats_snippets] , ecolor='k', color ='k',ls='none')
    ax.errorbar(Ra_eff_snippets,[snip[0] for snip in conduit_stats_snippets],yerr=[snip[1] for snip in conduit_stats_snippets], ecolor='g', color='g',ls='none')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.grid()
    #plt.show()
    plt.clf()
    plt.close()

    RA_EFF_HEAD.append(Ra_eff_head_snippets)
    RA_EFF_COND.append(Ra_eff_cond_snippets)

    RA_EFF.append(Ra_eff_snippets)
    HEAD_STATS.append(head_stats_snippets)
    CONDUIT_STATS.append(conduit_stats_snippets)
    ALL_HEAD.append(all_heads)
    ALL_CONDUITS.append(all_conduits)












# ------------ ALL MODELS, ALL SNIPPETS PLOTS ---------------------------



ALL_HEAD_FLAT = list()
ALL_CONDUITS_FLAT =list()
ALL_RA_HEAD_FLAT = list()
ALL_RA_COND_FLAT = list()

fig, ax = plt.subplots(1)
fig.set_size_inches(5,3)

for mod_idx in range(len(RA_EFF)):

    ALL_HEAD_flat = [item for sublist in ALL_HEAD[mod_idx] for item in sublist]
    ALL_CONDUITS_flat = [item for sublist in ALL_CONDUITS[mod_idx] for item in sublist]
    ALL_RA_EFF_HEAD_flat = [item for sublist in RA_EFF_HEAD[mod_idx] for item in sublist]
    ALL_RA_EFF_COND_flat = [item for sublist in RA_EFF_COND[mod_idx] for item in sublist]
    ALL_HEAD_FLAT.extend(ALL_HEAD_flat)
    ALL_CONDUITS_FLAT.extend(ALL_CONDUITS_flat)
    ALL_RA_HEAD_FLAT.extend(ALL_RA_EFF_HEAD_flat)
    ALL_RA_COND_FLAT.extend(ALL_RA_EFF_COND_flat)


    Ra_eff_snippets = RA_EFF[mod_idx]
    head_stats_snippets = HEAD_STATS[mod_idx]
    conduit_stats_snippets = CONDUIT_STATS[mod_idx]
    all_heads = ALL_HEAD[mod_idx]
    all_conduits = ALL_CONDUITS[mod_idx]

    ax.scatter(Ra_eff_snippets,[snip[0] for snip in head_stats_snippets],color='k',s=100, marker='o')
    ax.scatter(Ra_eff_snippets,[snip[0] for snip in conduit_stats_snippets], color='g',s=100, marker ='o')
    
    # for sn_idx in range(len(all_heads)):
    #     ax.scatter([Ra_eff_snippets[sn_idx]]*len(all_heads[sn_idx]), all_heads[sn_idx], color ='k', alpha = 0.05)
    #     ax.scatter([Ra_eff_snippets[sn_idx]]*len(all_conduits[sn_idx]), all_conduits[sn_idx], color ='g', alpha = 0.05)
    
    ax.errorbar(Ra_eff_snippets,[snip[0] for snip in head_stats_snippets],yerr=[snip[1] for snip in head_stats_snippets] , ecolor='k', color ='k',ls='none')
    ax.errorbar(Ra_eff_snippets,[snip[0] for snip in conduit_stats_snippets],yerr=[snip[1] for snip in conduit_stats_snippets], ecolor='g', color='g',ls='none')

print(len(ALL_HEAD_FLAT), len(ALL_CONDUITS_FLAT), len(ALL_RA_HEAD_FLAT),len(ALL_RA_COND_FLAT))

b_head,a_head,rsqu_head, V_head = this_model.getPowerLawCoeffs (ALL_RA_HEAD_FLAT , ALL_HEAD_FLAT)
b_cond,a_cond,rsqu_cond, V_cond = this_model.getPowerLawCoeffs (ALL_RA_COND_FLAT , ALL_CONDUITS_FLAT)
powerline_head= np.e**(a_head) * np.array(ALL_RA_HEAD_FLAT)**b_head
powerline_cond= np.e**(a_cond) * np.array(ALL_RA_COND_FLAT)**b_cond

ax.plot(ALL_RA_HEAD_FLAT, powerline_head, color='gray')
ax.plot(ALL_RA_COND_FLAT, powerline_cond, color='gray')

print("---- w = A Ra ^ b -------")
print("A_head = ", np.e**(a_head), " +/-", np.e**(a_head) * np.sqrt(V_head[1][1]))
print("b_head = ", b_head, "+/-", np.sqrt(V_head[0][0]))
print("r^2 = ", rsqu_head)

print("A_cond = ", np.e**(a_cond), " +/-", np.e**(a_cond) * np.sqrt(V_cond[1][1]))
print("b_cond = ", b_cond, "+/-", np.sqrt(V_cond[0][0]))
print("r^2 = ", rsqu_cond)


# ------------- Plot lithgow-bertelloni -----------------

L_analog = 0.3 #m

# cm to m
w_head_LB = 2 * 8.6e-2 * (this_model.D/L_analog) * np.array(ALL_RA_HEAD_FLAT)**(-0.14)
w_cond_LB = 2 * 6.8e-2 * (this_model.D/L_analog) * np.array(ALL_RA_COND_FLAT)**(-0.21)

ax.plot(ALL_RA_HEAD_FLAT, w_head_LB, color='b')
ax.plot(ALL_RA_COND_FLAT, w_cond_LB, color='b',linestyle='--')

ax.set_yscale('log')
ax.set_xscale('log')
ax.grid()
plt.show()
plt.clf()
plt.close()








