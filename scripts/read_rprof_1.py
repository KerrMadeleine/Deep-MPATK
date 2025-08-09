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
from modelparameters import TSURF, CORE_RADIUS, DOMAIN_THICKNESS

rprof_n0 = 0 
rprof_n1 = 57

for m in range(len(STAGYY_MODS)):

    # if '1300' not in STAGYY_MODS[m]:
    #     continue

    path = STAGYY_OUTPUT_FOLDER + STAGYY_MODS[m]

    # define model
    test_model = smc.StagYYModel(path, MESH_DIR, STAGYY_MODS_YDIM[m], STAGYY_MODS_ZDIM[m], STAGYY_MODS_NMAX[m], STAGYY_MODS_TCMB[m], STAGYY_MODS_TM[m], TSURF, CORE_RADIUS, DOMAIN_THICKNESS)
    test_model.set_rprof_dictionary(rprof_n0,rprof_n1)
    data_dict = test_model.rprof_data_dict

    time = data_dict['time'][:,0]
    ntime = len(time)

    # color map plot of temperature deviation from the initial mean value
    fig, ax = plt.subplots(1)
    fig.set_size_inches(9,3)
    c=ax.pcolor(data_dict['time'],data_dict['r'],data_dict['Tmean']-test_model.Tm_0,cmap=cm.vik, vmin=-200,vmax=200) #define color bar bounds with vmin and vmax
    cax = plt.axes([0.92, 0.25, 0.02, 0.5])
    fig.colorbar(c,cax=cax)
    #plt.show()
    plt.clf()
    plt.close()

    # get Vrms in different regions
    wholedomain_vrms_LM = test_model.getVrmsevolution(2,6) # domain divided into eights
    wholedomain_vrms_WM = test_model.getVrmsevolution(0,8)
    wholedomain_vrms_UM = test_model.getVrmsevolution(7,8)


    # plot Vrms in different regions
    fig, ax = plt.subplots(1)
    fig.set_size_inches(9,3)
    ax.plot(time, wholedomain_vrms_LM,label='LM')
    ax.plot(time, wholedomain_vrms_WM,label='WM')
    ax.plot(time, wholedomain_vrms_UM,label='UM')
    plt.grid()
    ax.set_ylabel('vel (mm/yr)')
    ax.set_xlabel('time')
    #plt.show()
    plt.clf()
    plt.close()

    # print out the average vrms
    print(STAGYY_MODS[m], 'Vrms (final half):', np.mean(wholedomain_vrms_LM[ntime//2:]))
