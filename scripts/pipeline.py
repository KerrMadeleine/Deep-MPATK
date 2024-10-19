import set_DeepMPATK_pars as par
from matplotlib.cbook import flatten
from sklearn.cluster import KMeans
from sklearn import metrics as skmet;
import numpy as np
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import csv
import os
import sys
import stagyypythonmodule as spm
from cmcrameri import cm

'''
This script is run 3 times throughout the analysis process and a series of binary flags from 
the shell script indicate what inputs/outputs this pipeline produces. The first one is to read
in the models that have a high time resolution to get the point of convective instability at the 
bottom TBL.
The second run is to obtain information about the depth of the convective domain in the stagnant lid 
planet for the full models and write that information to file for easier use later.
The third run is to use that information to detect plume material.
'''

# TODO: move this to the end of the global parameters file
############### COLOR MAPS 
c_white = clr.colorConverter.to_rgba('white',alpha = 0)
c_red=clr.colorConverter.to_rgba('red',alpha = 1)
c_magenta= clr.colorConverter.to_rgba('magenta',alpha = 1)
c_orange= clr.colorConverter.to_rgba('orange',alpha = 1)
c_yellow=clr.colorConverter.to_rgba('yellow',alpha = 1)
c_yegr=clr.colorConverter.to_rgba('greenyellow',alpha = 1)
c_green=clr.colorConverter.to_rgba('green',alpha = 1)
c_blue=clr.colorConverter.to_rgba('blue',alpha = 1)
c_brown=clr.colorConverter.to_rgba('brown',alpha = 1)
c_purple=clr.colorConverter.to_rgba('purple',alpha = 1)
c_black=clr.colorConverter.to_rgba('black',alpha = 1)
cmap_rb = clr.LinearSegmentedColormap.from_list('rb_cmap',[c_white,c_red,c_orange,c_yellow,c_yegr,c_green,c_blue,c_purple,c_magenta,c_brown,c_black],11)
cmap_bin= clr.LinearSegmentedColormap.from_list('bin_color',[c_white,c_black],2)

###############

################################
########  DEFINE MODULES ###### 
################################
##################################################################################

DIR = par.STAGYY_OUTPUT_FOLDER
TIMEDIR = par.STAGYY_OUTPUT_FOLDER
MESHDIR = par.STAGYY_OUTPUT_MESH_FOLDER

MODS=par.STAGYY_OUTPUT_MODS
TCMB=par.STAGYY_MODS_TCMB
TM = par.STAGYY_MODS_TM
SPT0 = par.STAGYY_MODS_SPT0
SPT1 = par.STAGYY_MODS_SPT1
DIMS = par.STAGYY_MODS_DIMS

haveBLX_RAC=par.PIPELINE_HAS_BLX_RA

if haveBLX_RAC:
    BL_X_KM = par.STAGYY_MODS_BL_X_KM
    BL_X_PT = par.STAGYY_MODS_BL_X_PT
    RAC = par.STAGYY_MODS_RAC
    N_ONSET=par.STAGYY_MODS_N_ONSET
    TAU = par.STAGYY_MODS_TAU

KMAXNUM=par.PIPELINE_KMAXNUM
wallsizefactor=par.PIPELINE_wallsizefactor
calcKmeans=par.PIPELINE_calcKmeansBool
calcConvReg=par.PIPELINE_CALC_CONV_REG
BLplots=par.PIPELINE_makeBLplots

# parameter flags controlling visualizations and statistics
computecoldTBL=par.PIPELINE_computecoldTBL
makeRavTimePlot=par.PIPLELINE_makeRavTimePlot
makeBLvsTimeplot=par.PIPELINE_makeBLvsTimeplot
makeannotatedTprof=par.PIPELINE_makeannotatedTprof
plotallmodsRavTime=par.PIPELINE_plotallmodsRavTime
showStatsandScores=par.PIPELINE_showStatsandScores
showViscprof=par.PIPELINE_showViscprof

# parameter flags to indicate which path down the pipeline we traverse
Write_Conv_Reg_Thick=par.PIPELINE_WRITE_CONV_REG
Write_Bl_Rac_Pars=par.STAGYY_WRITE_BL_RAC_PARS
Write_Raeff=par.STAGYY_WRITE_RA_EFF

print(Write_Conv_Reg_Thick, Write_Bl_Rac_Pars, Write_Raeff)

startindmod=par.STAGYY_OUTPUT_MODS_0
endin=par.STAGYY_OUTPUT_MODS_1

##############################################################################

MODS=MODS[startindmod:endin]
TCMB=TCMB[startindmod:endin]
TM=TM[startindmod:endin]
SPT0=SPT0[startindmod:endin]
SPT1=SPT1[startindmod:endin]
DIMS=DIMS[startindmod:endin]

if haveBLX_RAC:
    RAC=RAC[startindmod:endin]
    BL_X_KM=BL_X_KM[startindmod:endin]
    BL_X_PT=BL_X_PT[startindmod:endin]
    N_ONSET=N_ONSET[startindmod:endin]
    TAU=TAU[startindmod:endin]
else:
    TAU=np.ones(len(MODS))

print(MODS) # PRINT STATEMENT

DIRS=[DIR+mod for mod in MODS]
TDIRS=[TIMEDIR+mod for mod in MODS]
# NUM_TimeSteps (the exact number of output files output by StagYY)
NUM_TimeSteps = [len([file for file in os.listdir(dir) if file[0:2]=='t_']) for dir in DIRS] #new
# Turn the number of output files into a string to loop through

NUMS=[[str(i).zfill(5) for i in range(0,NUM_TimeSteps[j])] for j in range(len(MODS))]

##############################################################################
# GLOBAL VARS
##############################################################################
domain_THICK = 2942.e3
air_THICK = 0
core_THICK = 3110.e3

E_act=300e3
A_arrhenius=1e20/np.exp(E_act/(8.314*1600))

a_exp = 3.e-5 # K^-1
g=8.87 # m/s^2
Cp=1200 
kappa=1.e-6 #m^2/s
rho=3300
therm_conductivity = kappa*rho*Cp
H_flux = 0 #1.623e-8 # W/m^3
H_mass = 0 #4.9185e-12 # W/kg

D_MantlenStag = domain_THICK-air_THICK
mantle_THICK=domain_THICK-air_THICK
ra_num_const = (mantle_THICK)**3 * rho**2 * Cp*g*a_exp/therm_conductivity
mantle_THICK=mantle_THICK/1000

geom_pars = [domain_THICK,air_THICK,core_THICK,D_MantlenStag,mantle_THICK]
phys_pars=[a_exp,g,Cp,kappa,rho,therm_conductivity,H_flux,H_mass,ra_num_const]

##############################################################################



##############################################################################
# LOOP THROUGH DIRECTORIES
##############################################################################

BLTHICKkm_itpl=[]
BLTHICKp_itpl=[]

TIMES=[]
ALL_VMEAN=[]
ALL_TMEAN=[]
ALL_DELT = []

T_for_kellog = []

RA=[]
RA_BIN=[]

PLMNO=[]
PLMTK=[]
NO_TK_plumes=[]

CONVREGTHK_ALL=[]
BOTCLDTBLTEMP_ALL=[]
FULLCONV_RABOT_ALL=[]
FULLCONV_RAINT_ALL=[]
FULLCONV_RACOMB_ALL=[]
TIMESCONVRA_ALL=[]

if Write_Raeff:
    raeff_file = open(par.RA_EFF_FILE+'raeff.txt', 'w')
    raeff_file.close()
    raeff_file = open(par.RA_EFF_FILE+'raeff.txt', 'a')

for i in range(len(DIRS)): #for each of the models

    dir=DIRS[i]
    nums=NUMS[i]
    t0 = TM[i]
    tcmb=TCMB[i]
    spt0=SPT0[i]
    spt1=SPT1[i]

    Times=spm.getTimeARR(TDIRS[i],spt0,spt1,1)  #time directory, step size

    if Write_Bl_Rac_Pars:

        onset_file = open(par.CONV_ONSET_FILE+MODS[i][:-1]+'.txt', 'w')
        onset_file.close()
        onset_file = open(par.CONV_ONSET_FILE+MODS[i][:-1]+'.txt', 'a')

    if haveBLX_RAC:

        n_onset=N_ONSET[i]
        print('n onset:',n_onset)
        spt0=int(np.ceil(n_onset))
        Times=spm.getTimeARR(TDIRS[i],spt0,spt1,1)
        Ra_c=RAC[i]
        blX_km=BL_X_KM[i]
        blX_pt=BL_X_PT[i]


    meshdim=DIMS[i]
    if meshdim=='1024x256':
        mesh_y_dim=1024
    elif meshdim=='512x128':
        mesh_y_dim=512
    elif meshdim=='256x64':
        mesh_y_dim=256
    elif meshdim=='400x100':
        mesh_y_dim=400
    else:
        print('MESH ERROR')
        break

    ######################################################################################
    ##### MESHES
    ################################################################################
    #OPEN MESH FILE	
    Tablemesh=[]

    with open('/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/STAGmesh'+meshdim+'.csv') as f:
        reader = csv.reader(f)
        for col in reader:
            Tablemesh.append(col)
                    
    y=[float(i) for i in Tablemesh[0][0:mesh_y_dim]] # get the radii coords
    z=[float(i) for i in Tablemesh[0][mesh_y_dim::2*mesh_y_dim]] # get the z coords

    Y=len(y)
    Z=len(z)

    arrY=[i for i in range(1,Y+1)]
    arrZ=[i for i in range(1,Z+1)]

    eighth_of_dom=(mesh_y_dim//4)//8
    mm0=2*eighth_of_dom
    mm1=5*eighth_of_dom

    if haveBLX_RAC:
        WALLSIZE=int(np.floor(wallsizefactor*blX_pt)) #mesh_y_dim//128
    else:
        WALLSIZE=0

    WALLBOUNDL=WALLSIZE
    WALLBOUNDR=Y-WALLSIZE

    midmantlen = [mm0,mm1]
    dims = [Y,Z]
    # (DEFINED ABOVE) geom_pars = [domain_THICK,air_THICK,core_THICK,D_MantlenStag,mantle_THICK]
    # (DEFINED ABOVE) phys_pars=[a_exp,g,Cp,kappa,rho,therm_conductivity,H_flux,H_mass,ra_num_const]

    pars = [dims,phys_pars,geom_pars,midmantlen]
    RAD = spm.getRAD(TIMEDIR+MODS[i],pars)

    print('******** CURRENT DIR:',dir,' *********') # PRINT STATEMENT

    BLTHICKkm_1mod_itpl=[]
    BLTHICK_follow = []
    BLTHICKp_1mod_itpl=[]
    RA_BL_1mod=[]
    BL3_eta_1mod=[]
    vmean_1mod=[]
    delt_1mod=[]

    vmean_truevtrim_1mod=[]
    tmean_truevtrim_1mod=[]

    TIMES.append(Times)
    TPlot_for_kellog=[]
    RAplot=[]
    
    RAcritBIN=[]
    RAgrpPLT=[]

    TMEAN_1mod=[]
    PLMNO_1mod=[]
    PLMTK_1mod=[]
    NO_TK_plumes_1mod=[]

    # initialize BL thickness, Ra_crit, and t_onset variables with min values
    BLDIM_p=0
    BLDIM_km=0.1
    ETA_half=1e25
    DELT_rac=1
    RA_CRIT=1000
    t_onset=0

    if computecoldTBL:
        COLDBLTHICK_1mod=[]
        TPlot_for_kellog_UP_TBL=[]

    CONVREGTHK=[]
    BOTCLDTBLTEMP=[]
    FULLCONV_RABOT=[]
    FULLCONV_RAINT=[]
    FULLCONV_RACOMB=[]
    TIMESCONVRA = []
    TIMESTAURA =[]

    #################################
    for n in nums[spt0:spt1]: #for each time-step of each model

        print('TIMESTEP:  ',n) # PRINT STATEMENT
        N=int(n)
        #writerAVG.writerow([t]) # add the model time to the output csv file

        # READ THE MOD/TS FILE
        openFileT= dir+'t_'+n+'.csv' #temp file to open
        openFileV= dir+'eta_'+n+'.csv' #visc file to open

        #OPEN TEMP FILE
        TableT=[]
        with open(openFileT) as f:
            reader = csv.reader(f,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)
            for col in reader:
                #print(col)
                TableT.append(col)
                
        #OPEN VISC FILE	
        TableV=[]	
        with open(openFileV) as f:
            reader = csv.reader(f, delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)
            for col in reader:
                #print(col)
                TableV.append(col)


        TableT=np.ravel(TableT,order='C')
        TableV=np.ravel(TableV,order='C')

        TGRAD=[spm.trimmedMean(spm.getRadSlice(i,TableT,pars)) for i in range(Z)]
        VGRAD=[spm.trimmedMeanV(spm.getRadSlice(i,TableV,pars)) for i in range(Z)]

        TGRAD_TRUE=[np.mean(spm.getRadSlice(i,TableT,pars)) for i in range(Z)]
        VGRAD_TRUE=[np.exp(np.mean(np.log(spm.getRadSlice(i,TableV,pars)))) for i in range(Z)]

        TGRAD=np.array(TGRAD)
        VGRAD=np.array(VGRAD)

        TMEAN = np.mean(TGRAD[mm0:mm1])
        VMEAN = np.mean(VGRAD[mm0:mm1])

        vmean_1mod.append(VMEAN)
        TMEAN_1mod.append(TMEAN)
        tmean_truevtrim_1mod.append(TGRAD_TRUE-TGRAD)
        vmean_truevtrim_1mod.append(VGRAD_TRUE-VGRAD)

        minVISC=np.min(TableV) 

        #print('RATIO (vmean, vmin, vmin, rat1, rat2):', VMEAN, np.min(VGRAD), minVISC, VMEAN/np.min(VGRAD),VMEAN/minVISC)

        if calcConvReg:

            convregthk, botstag_temp = spm.getSTAGTHICK(TGRAD,VGRAD,RAD,Times[N-spt0],N,MODS[i][:-1],pars,tcmb,'trim',showViscprof,par.PIPELINE_VprofFreq)
            delTeff = rho*H_mass*(convregthk*10**3)**2/therm_conductivity

            rabot = (convregthk*10**3)**3 * rho**2 * Cp * g * a_exp * (tcmb-botstag_temp) / (therm_conductivity * VMEAN)
            raint = (convregthk*10**3)**3 * rho**2 * Cp * g * a_exp * delTeff / (therm_conductivity * VMEAN)
            racomb = (convregthk*10**3)**3 * rho**2 * Cp * g * a_exp * (delTeff+(tcmb-botstag_temp))/ (therm_conductivity*VMEAN)
            
            CONVREGTHK.append(convregthk)
            BOTCLDTBLTEMP.append(botstag_temp)

            FULLCONV_RABOT.append(rabot)
            FULLCONV_RAINT.append(raint)
            FULLCONV_RACOMB.append(racomb)

            TIMESCONVRA.append(Times[N-spt0])
            TIMESTAURA.append(Times[N-spt0]/TAU[i])

            #print('ts:',N,' Ra_bot:',rabot,' Ra_int',raint,' Ra_combined:',racomb)

        delt=tcmb-TMEAN
        delt_1mod.append(delt)

        ########################################

        THIGH=TMEAN+(0.1*delt)                                  # % TBL TEMP VALUE

        if computecoldTBL and calcConvReg:

            TLOW=TMEAN-(0.1*(TMEAN-botstag_temp)) #prof,rad,t,tstag,convregthk

            [BLkm_itpl_cold,BLp_itpl_cold] = spm.getBLptsITPL_cold(TGRAD,RAD,TLOW,botstag_temp,convregthk,tcmb)
        ### FROM TURCOTTE AND SHUBERT PG 184 ON HEAT TRANSFER
# ******************************************************************************************************************************************************
        # BL across Y
        BLS_kms=[spm.getBLptsITPL(spm.getTprofile(yind,TableT,pars),RAD,THIGH,tcmb)[0] for yind in range(Y)]
        BLS_pts=[spm.getBLptsITPL(spm.getTprofile(yind,TableT,pars),RAD,THIGH,tcmb)[1] for yind in range(Y)]
        
        #POLY FIT across Y
        poly_km = np.polyfit(arrY, BLS_kms, deg=2)
        poly_pts = np.polyfit(arrY, BLS_pts, deg=2)        
        BL_poly_km=np.polyval(poly_km, arrY)
        BL_poly_pt=np.polyval(poly_pts, arrY)
        
        # single value BL calc from TGRAD
        [BLkm_itpl,BLp_itpl] = spm.getBLptsITPL(TGRAD,RAD,THIGH,tcmb)

        # the final value for this timestep:
        if BLDIM_km<BLkm_itpl:

            BLDIM_km=BLkm_itpl
            BLDIM_p=BLp_itpl
            ETA_half=10**(spm.getInterpVal(BLDIM_km,np.log10(VGRAD),RAD,np.log10(A_arrhenius*np.exp(E_act/(8.314*tcmb)))))
            DELT_rac=delt
            RA_CRIT=(a_exp*rho*g*DELT_rac*(BLDIM_km*10**3)**3/(ETA_half*kappa))
            ts_onset=N
            t_onset=Times[N-spt0]


        eta_half=10**(spm.getInterpVal(BLkm_itpl,np.log10(VGRAD),RAD,np.log10(A_arrhenius*np.exp(E_act/(8.314*tcmb)))))

        #OVERWRITING BL with EXPIRICAL VALUES FROM OTHER EXP.
        if haveBLX_RAC:

            BLDIM_km=blX_km
            BLDIM_p=blX_pt
            RA_CRIT=Ra_c

        BLTHICKkm_1mod_itpl.append(BLDIM_km)
        BLTHICK_follow.append(BLkm_itpl)

        if haveBLX_RAC:
            RA_BL_1mod.append(Ra_c)
        else:
            RA_BL_1mod.append(a_exp*rho*g*delt*(BLkm_itpl*10**3)**3/(eta_half*kappa))

        BL3_eta_1mod.append((BLkm_itpl*1000)**3/VMEAN)

        Tshell_for_kellog = spm.getRadSlice_ITPL(BLDIM_p,0,TableT,pars)
        Vshell_for_kellog = spm.getRadSlice_ITPL(BLDIM_p,0,TableV,pars)

        RAloc_2 = spm.getRAlslice(BLDIM_p,BLDIM_km,delt,TableV,pars) #polynomial
        RAloc_2=RAloc_2[WALLBOUNDL:WALLBOUNDR]

        RAloc = spm.getRAlslice(BLp_itpl,BLkm_itpl,delt,TableV,pars)
        RAloc=RAloc[WALLBOUNDL:WALLBOUNDR]
        RAplot.append(RAloc_2)

        if computecoldTBL and haveBLX_RAC:

            print('time:',Times[N],' \n tau:', Times[N]/TAU[i], '\nTMEAN:', TMEAN, '\nT_l',THIGH, '\ntcmb',tcmb,'\ndelT:', delt,'\ntlith',botstag_temp, '\neta:', VMEAN,'\nBL:',BLkm_itpl)
        
        if Write_Bl_Rac_Pars:

            if (N-spt0)>=1:
                if BLTHICKkm_1mod_itpl[N]!=BLTHICKkm_1mod_itpl[N-1]:
                    onset_file.write('\nBL, '+ str(BLDIM_km) + ',' +str(BLDIM_p) + ','+ str(RA_CRIT))

            else:
                onset_file.write('\nBL, '+ str(BLDIM_km) + ',' +str(BLDIM_p) + ','+ str(RA_CRIT))
        
        if BLplots and makeBLvsTimeplot and N%par.PIPELINE_BLplotfreq==0:

            fig, ax = plt.subplots(1)
            fig.set_size_inches(7,5)
            fig.suptitle(MODS[i][:-1])

            ax.plot(np.array(Times[:len(BLTHICKkm_1mod_itpl)])/TAU[i],BLTHICK_follow)
            ax.plot(np.array(Times[:len(BLTHICKkm_1mod_itpl)])/TAU[i],BLDIM_km*np.ones(len(BLTHICKkm_1mod_itpl)))
            
            ax.set_ylabel('$l$ (km)',fontsize=15)
            ax.set_xlabel('overturns',fontsize=15)
            ax.grid(axis='both', which='major', color='0.25', linestyle='-')
            ax.grid(axis='x', which='minor', color='0.8', linestyle='--')
            plt.show()
            #plt.savefig('/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/PLOTS/JUN22/WALL_TBLthk/TBL_crit/'+MODS[i][:-1]+'_'+str(N)+'.png',dpi=300)
            plt.clf()
            plt.close()


            
            fig, ax = plt.subplots(1)
            fig.set_size_inches(10,7)
            fig.suptitle(MODS[i][:-1])

            ax.plot(np.arange(len(RAloc_2)),np.array(RAloc_2),color='k')
            ax.plot(np.arange(len(RAloc)),np.array(RAloc),color='blue')
            ax.plot(np.arange(Y),np.ones(Y)*RA_CRIT,color='red')
            
            ax.set_xticks(ticks=[0,Y//2,Y],labels=[0,'$\pi/2$','$\pi$'],fontsize=15)
            ax.set_yscale('log')
            ax.set_ylabel('RaC and RaL')
            ax.grid(axis='both', which='major', color='0.25', linestyle='-')
            ax.grid(axis='both', which='minor', color='0.8', linestyle='--')
            plt.show()
            #plt.savefig('/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/PLOTS/RA_C/'+MODS[i][:-1]+'_'+str(N)+'.png',dpi=300)
            plt.clf()
            plt.close()

        TPlot_for_kellog.append(Tshell_for_kellog[WALLBOUNDL:WALLBOUNDR])
        ######################################

        if makeannotatedTprof and N%par.PIPELINE_TprofFreq==0: # plot the boundary layer thickness across the Y domain
            
            fig, ax = plt.subplots(1)
            fig.set_size_inches(6,5)
            fig.suptitle(MODS[i][:-1]+'bl:'+str(np.mean(BLDIM_km))+'time:'+str(Times[N-spt0]/TAU[i]))
            ax.plot(TGRAD,RAD,color='k',linewidth=5)

            ax.plot(np.ones(8)*THIGH,np.linspace(0,2000,8),color='green')
            ax.plot(arrY,[BLkm_itpl for ii in range(Y)])
            ax.plot([600,2100],[BLkm_itpl,BLkm_itpl],color='green',linewidth=2,linestyle='dashed')
            ax.plot([THIGH,THIGH],[0,3000],color='green',linewidth=2)
            ax.plot([tcmb,tcmb],[0,3000],color='red',linewidth=2)
            if par.PIPELINE_computecoldTBL:
                ax.plot([botstag_temp,botstag_temp],[0,3000],color='purple',linewidth=2)

            ax.set_yticks(ticks=np.arange(0,3000,500),labels=np.arange(0,3000,500),fontsize=15)
            ax.set_ylim([0,mantle_THICK])
            ax.set_xlim([700,tcmb+100])
            ax.set_xticks(ticks=np.arange(700,tcmb+100,200),labels=np.arange(700,tcmb+100,200),fontsize=15)
            ax.grid(axis='both', which='major', color='0.25', linestyle='-')
            ax.grid(axis='x', which='minor', color='0.8', linestyle='--')
            #plt.savefig('/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/PUB_FIG/t_prof/'+MODS[i]+str(N)+'.png',dpi=300)
            plt.show()
            plt.clf()
            plt.close()

        if haveBLX_RAC==0:
            continue

        Ra_crit=RA_CRIT
        RA_binSLC=[el>=Ra_crit for el in RAloc_2]

        RAcritBIN.append(RA_binSLC)  # change here for plume vs slab detection
# ******************************************************************************************************************************************************

        [pnum,pwid]=spm.plumeNo_Th(Ra_crit,BLkm_itpl,RAloc_2,1,pars)
        #[pnum,pwid]=spm.plumeNo_Th(Ra_crit,BLkm_itpl,RAloc,1,pars) #REMOVE THIS

        NO_TK_plumes_1mod.append([pnum,pwid])
        #trim_wid_lst=[wid for wid in pwid if (wid<=1000 and wid>40)]
        trim_wid_lst=[wid for wid in pwid]     

        if pnum==0:
            PLMNO_1mod.append(0)
            PLMTK_1mod.append([])
        else:
            PLMTK_1mod.append(trim_wid_lst)
            PLMNO_1mod.append(len(trim_wid_lst)) 
# ******************************************************************************************************************************************************
    if calcConvReg and makeRavTimePlot:

        fig, ax = plt.subplots(3)
        #print([np.mean(FULLCONV_RABOT[i:i+20]) for i in range(0,len(CONVREGTHK)-21,20)])
        ax[0].scatter(TIMESTAURA,CONVREGTHK)
        ax[1].scatter(TIMESTAURA,[tcmb-el for el in BOTCLDTBLTEMP])
        ax[2].scatter(TIMESTAURA,FULLCONV_RABOT)
        ax[2].scatter(TIMESTAURA,FULLCONV_RAINT)
        ax[2].scatter(TIMESTAURA,FULLCONV_RACOMB)

        ax[2].set_xlabel('Time (overturns)')
        ax[0].set_ylabel('Conv. Dom. (km)')
        ax[1].set_ylabel('$\Delta T$ dom.')
        ax[2].set_ylabel('Ra')
        ax[2].set_yscale('log')
        print(np.mean(FULLCONV_RABOT))
        print(np.std(FULLCONV_RABOT))

        plt.show()
        plt.clf()
        plt.close()
    
    if Write_Conv_Reg_Thick:

        onset_file = open(par.CONV_ONSET_FILE+MODS[i][:-1]+'.txt', 'a')
        onset_file.write('\n '+str(len([BLTHICKkm_1mod_itpl[nn] for nn in range(len(BLTHICKkm_1mod_itpl)-1) if (BLTHICKkm_1mod_itpl[nn]<BLTHICKkm_1mod_itpl[nn+1])])))
        onset_file.write('\n ConvRegTh:, '+ str(np.mean(CONVREGTHK[30:])))
        onset_file.close()

    ####################################
    print('TBL THICKNESS:',BLDIM_km) # PRINT STATEMENT

    BLTHICKkm_itpl.append(BLTHICKkm_1mod_itpl)
    BLTHICKp_itpl.append(BLTHICKp_1mod_itpl)

    ALL_VMEAN.append(vmean_1mod)
    ALL_DELT.append(delt_1mod)
    ALL_TMEAN.append(TMEAN_1mod)
    NO_TK_plumes.append(NO_TK_plumes_1mod)
    PLMNO.append(PLMNO_1mod)

    allplumes_1mod=sorted(flatten(PLMTK_1mod))
    allplumes=np.array(allplumes_1mod)
    allplumes.reshape(-1,1)
    
    PLMTK.append(sorted(allplumes))

    CONVREGTHK_ALL.append(CONVREGTHK)
    BOTCLDTBLTEMP_ALL.append(BOTCLDTBLTEMP)
    FULLCONV_RABOT_ALL.append(FULLCONV_RABOT)
    FULLCONV_RAINT_ALL.append(FULLCONV_RAINT)
    FULLCONV_RACOMB_ALL.append(FULLCONV_RACOMB)
    TIMESCONVRA_ALL.append(TIMESCONVRA)

    if Write_Raeff:
        raeff_file.write(str(FULLCONV_RABOT_ALL))

    if Write_Bl_Rac_Pars:
        onset_file.close()
    if calcKmeans:

        #  * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        # * * * * * * * *     KMEANS      * * * * * * * * * * * * 
        # * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        kmax=KMAXNUM
        CENTROIDS_1mod=[]
        print('NUMBER OF PLUMES:',len(allplumes))
        data = allplumes.reshape(-1,1)
        #data = [radius**2 for radius in data]
        stats_empty=1

        if len(allplumes)>=kmax:

            SSdiff_CC=[]
            silo_1mod=[]
            DB_1mod=[]
            CH_1mod=[]
            GAP_1mod=[]
            gapstat_1mod=[]
            SSdiff_elbow=[]

            for ii in range(1,kmax):

                kmeans = KMeans(n_clusters=ii, random_state=0).fit(data)

                ssdiff=kmeans.inertia_
                labels=kmeans.labels_
                std_per_feat=[]
                num_per_feat=[]

                for jj in range(ii):
                    std_per_feat.append(np.std([allplumes[x] for x in range(len(allplumes)) if labels[x]==jj]))
                    num_per_feat.append(len([allplumes[x] for x in range(len(allplumes)) if labels[x]==jj]))
                
                SSdiff_CC.append(ssdiff)
                SSdiff_elbow.append([ii,ssdiff])

                if ii!=1:

                    silo=skmet.silhouette_score(data, labels, metric='euclidean')
                    silo_1mod.append(silo)

                    DBsc=skmet.davies_bouldin_score(data, labels)
                    DB_1mod.append(DBsc)
                    stats_empty=0

                    CHsc=skmet.calinski_harabasz_score(data, labels)
                    CH_1mod.append(CHsc)

                else:
                    silo_1mod.append(0)
                    DB_1mod.append(1)
                    CH_1mod.append(0)
                    stats_empty=0

                refDisps=[]
                for r in range(100):

                    # Create new random reference set
                    randomReference = np.random.random_sample(size=data.shape)
                    # Fit to it
                    km = KMeans(ii)
                    km.fit(randomReference)
                    refDisp = km.inertia_
                    refDisps.append(refDisp)
                
                # Fit cluster to original data and create dispersion
                origDisp = kmeans.inertia_
                
                # Calculate gap statistic
                gap = np.log(np.mean(refDisps)) - np.log(origDisp)
                # Assign this loop's gap statistic to gaps
                GAP_1mod.append(gap)

                means = kmeans.cluster_centers_
                means=flatten(means)
                means=list(means)

                #print(ii,":  ",means, std_per_feat, num_per_feat)

                clumpdata=list(zip(means,std_per_feat,num_per_feat))
                clumpdata=sorted(clumpdata)
                
                CENTROIDS_1mod.append([ii,[x[0] for x in clumpdata]])
                print('[',ii,',',[x[0] for x in clumpdata],',',[x[1] for x in clumpdata],',',[x[2] for x in clumpdata],']')
                
                if ii>1:
                    fig, ax = plt.subplots(1)
                    fig.set_size_inches(5,3)
                    fig.suptitle(MODS[i][:-1])
                    ax.scatter([x[0] for x in clumpdata],[x[2] for x in clumpdata])
                    ax.errorbar([x[0] for x in clumpdata],[x[2] for x in clumpdata],xerr=[x[1] for x in clumpdata])
                    #plt.savefig(str(ii)+'_'+MODS[i][:-1]+'line.png',dpi=300)
                    #plt.show()
                    plt.clf()
                    plt.close()

            SSdiff_elbow=np.array(SSdiff_elbow)
        else:
            print('No plumes')
            SSdiff_CC=np.arange(1,kmax,1)
        # ******************************************************************************************************************************************************

    T_for_kellog.append(TPlot_for_kellog)
    RA.append(RAplot)
    RA_BIN.append(RAcritBIN)#[spt0:spt1]) #SPLIT


# ******************************************************************************************************************************************************
    if calcKmeans:
        # USE RAcritBIN, NO_TK_plumes_1mod, CENTROIDS_1mod, low:30 and high:1000
        RAcritBIN=RAcritBIN #[spt0:spt1] # SPLIT
        testR=0
        arrY1=arrY
        arrY=np.array(arrY[WALLBOUNDL:WALLBOUNDR+1])*(np.pi/Y)
        if len(CENTROIDS_1mod)!=0:
            testR=1
            for cc in range(1,len(CENTROIDS_1mod),1):
                HEAD_COORDS=[]
                RAgrpPLT=[]
                print('GROUP K=',cc)
                for tt in range(len(RAcritBIN)): #JAN5
                    intMask=[int(x) for x in RAcritBIN[tt]]
                    newbinSLC,Y_HEAD_COORDS=spm.clumpedSlice(intMask,NO_TK_plumes_1mod[tt][0],NO_TK_plumes_1mod[tt][1],CENTROIDS_1mod[cc-1][0],CENTROIDS_1mod[cc-1][1])
                    RAgrpPLT.append(newbinSLC)
                    HEAD_COORDS.append([(np.pi*(Y_HEAD_COORDS[i]+WALLBOUNDL)/mesh_y_dim,Times[tt]) for i in range(len(Y_HEAD_COORDS))])
                RApts=RA[i] #[spt0:spt1] # SPLIT
                RAbinmask=RA_BIN[i]
                fileName=MODS[i][:-1]

                if par.PIPELINE_plotKmeansKellogg:
                    fig, ax = plt.subplots(1) # JAN4
                    fig.set_size_inches(16,7)
                    fig.suptitle(fileName+"  RA loc: k="+str(cc))
                    ax.set_ylabel('Time (Myrs)')
                    ax.set_xlabel('$\Phi$ (angles span core-mantle boundary)',size=20)
                    RApts=np.array(RApts)
                    RAgrpPLT=np.array(RAgrpPLT)
                    Times_yaxis=np.array(TIMES[i])
                    XCMB,TIMEevo = np.meshgrid(arrY,TIMES[i])
                    XCMB_grp,TIMEevo_grp = np.meshgrid(arrY,TIMES[i])
                    ax.pcolor(XCMB,TIMEevo,RApts,cmap='Blues',vmin=0, vmax=4000) #SPLIT
                    c=ax.pcolor(XCMB_grp,TIMEevo_grp,RAgrpPLT,cmap=cmap_rb,vmin=0,vmax=cc) # groups #SPLIT
                    ax.set_ylim([TIMES[i][0],TIMES[i][spt1-spt0]+1])
                    #ax.set_xlim(arrY1[0],arrY1[-1])
                    cax = plt.axes([0.93, 0.1, 0.02, 0.8])
                    #plt.savefig('/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/PLOTS/JUN22/WALL_TBLthk/K/'+MODS[i][:-1]+'/k'+str(cc)+'.png',dpi=300)
                    plt.show()
                    plt.clf()
                    plt.close()
        arrY=arrY1

    # ******************************************************************************************************************************************************

    ##############################################################################
    # PLOTS (each model)
    ##############################################################################
    if calcKmeans and showStatsandScores:
        fig, ax = plt.subplots(6)
        fig.set_size_inches(6,7)
        fig.suptitle(MODS[i][:-1])
        ax[0].plot(np.arange(1,kmax,1),SSdiff_CC)
        ax[1].plot(np.arange(1,kmax,1),SSdiff_CC)
        if stats_empty==0:
            ax[2].plot(np.arange(1,kmax,1),DB_1mod)
            ax[3].plot(np.arange(1,kmax,1),CH_1mod)
            ax[4].plot(np.arange(1,kmax,1),silo_1mod)
            ax[5].plot(np.arange(1,kmax,1),GAP_1mod)
        ax[1].set_yscale('log')
        ax[0].set_ylabel('ss diff')
        ax[1].set_ylabel('ss diff -log')
        ax[2].set_ylabel('DB')
        ax[3].set_ylabel('CH')
        ax[4].set_ylabel('silohuette')
        ax[5].set_ylabel('gapstat')
        ax[0].set_xlabel('k')
        plt.subplots_adjust(bottom=0.1, right=0.9, left=0.18, top=0.9, hspace=0.4)
        #plt.savefig('/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/PLOTS/JUN22/WALL_TBLthk/K/'+MODS[i][:-1]+'/kOPT.png',dpi=300)
        plt.show()
        plt.clf()
        plt.close()

    
        # * * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * *
    #  
    #  
    #                      TEMP/RA  KELLOGG
    TempKellogPlot=par.PIPELINE_tempKellogg
    if TempKellogPlot:
        print(len(T_for_kellog[i]),len(ALL_TMEAN[i]))
        Tpts=np.array(T_for_kellog[i]) #- np.array(ALL_TMEAN[i])[:, np.newaxis]
        if par.PIPELINE_tempKellogg_minusMean:
            Tpts=np.array(T_for_kellog[i]) - np.array(ALL_TMEAN[i])[:, np.newaxis]
        arrY1=np.insert(np.array(arrY[WALLBOUNDL:WALLBOUNDR])*(np.pi/Y), 0, 0)
        fileName=MODS[i][:-1]
        Times_yaxis=np.array(TIMES[i])

        fig, ax = plt.subplots(1)
        fig.set_size_inches(7,4)
        fig.suptitle('Temperature for model ' + fileName)
        print('arrY:',arrY1[0],arrY1[-1])
        print('times:',TIMES[i][0],TIMES[i][-1],spt0,spt1)
        XCMB,TIMEevo = np.meshgrid(arrY1,Times_yaxis/TAU[i])
        print(np.shape(XCMB),np.shape(TIMEevo),np.shape(Tpts[:-1]) )
        c=ax.pcolor(XCMB,TIMEevo,Tpts,cmap=cm.vik,vmin=-100,vmax=100)
        ax.set_xticks(ticks=[0,np.pi/4,np.pi/2,3*np.pi/4,np.pi],labels=[0,'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'],fontsize=15)
        if TAU[i]==1:
            ax.set_yticks(ticks=np.arange(0,2000,200),labels=np.arange(0,2000,200),fontsize=15)
        else:
            ax.set_yticks(ticks=np.arange(0,50,5),labels=np.arange(0,50,5),fontsize=15)
        ax.set_ylim([TIMES[i][0]/TAU[i],TIMES[i][-1]/TAU[i]])
        ax.set_xlabel('$\Phi$ (angles span core-mantle boundary)',size=20)
        ax.set_ylabel('# overturns',size=20)
        cax = plt.axes([0.93, 0.1, 0.02, 0.8])
        cax.set_label('$T$ (K)')
        fig.colorbar(c,cax=cax)
        #plt.savefig('/Users/mkerr/VenusResearch/Kellog-Clustering/scripts/APR7_all/temp/'+MODS[i][:-1]+'.png',dpi=300)
        plt.show()
        plt.clf()
        plt.close()

    #* * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * *


    # * * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * *
    #  
    #                       RAloc with Ra_crit mask
    testR=par.PIPELINE_plotRaloc_w_mask
    if testR==1:
        RApts=RA[i]
        RAbinmask=RAcritBIN
        arrY1=np.insert(np.array(arrY[WALLBOUNDL:WALLBOUNDR])*(np.pi/Y), 0, 0)
        fileName=MODS[i][:-1]
        Times_yaxis=np.array(TIMES[i])


        fig, ax = plt.subplots(1)
        fig.set_size_inches(15,7)
        fig.suptitle(fileName)
        XCMB,TIMEevo = np.meshgrid(arrY1,Times_yaxis/TAU[i])    # TIME- STEP
        c=ax.pcolor(XCMB,TIMEevo,RApts,cmap='Blues',vmin=0, vmax=4000)
        ax.plot(arrY1,5.447*np.ones(len(arrY1)),color='k')
        if TAU[i]==1:
            ax.set_yticks(ticks=np.arange(0,2000,200),labels=np.arange(0,2000,200),fontsize=15)
        else:
            ax.set_yticks(ticks=np.arange(0,50,5),labels=np.arange(0,50,5),fontsize=15)
        ax.set_xlabel('$\Phi$ (angles span core-mantle boundary)',size=20)
        #ax.set_ylabel('time (Myrs)',size=20)
        ax.set_ylabel('# overturns',size=20)
        ax.set_ylim([TIMES[i][0]/TAU[i],TIMES[i][-1]/TAU[i]])
        ax.set_xticks(ticks=[0,np.pi/4,np.pi/2,3*np.pi/4,np.pi],labels=[0,'$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'],fontsize=15)
        cax = plt.axes([0.93, 0.1, 0.02, 0.8])
        cax.set_label('$log10(Ra_l)$')
        fig.colorbar(c,cax=cax)
        #plt.savefig('/Users/mkerr/VenusResearch/Kellog-Clustering/scripts/APR7_all/kluster_regime4_smooth/binTS_'+MODS[i][:-1]+'.png',dpi=300)
        plt.show()
        plt.clf()
        plt.close()
    #* * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * ** * *

    ##############################################################################

if Write_Raeff:
    raeff_file.close()
##############################################################################
# PLOTS (all models)
##############################################################################

if calcConvReg and plotallmodsRavTime:

    fig, ax = plt.subplots(1)
    fig.set_size_inches(16,7)
    for j in range(len(MODS)):
        ax.plot(TIMESCONVRA_ALL[j],FULLCONV_RABOT_ALL[j])
    ax.set_ylabel('Ra bottom')
    ax.set_xlabel('time (Myr)')
    ax.set_yscale('log')
    #plt.savefig('/Users/mkerr/VenusResearch/2022/sept2022_analysis/OCT4_l-Ra-slab_vs_time/Ra_eff.png',dpi=400)
    plt.show()
    plt.clf()
    plt.close()

    fig, ax = plt.subplots(1)
    fig.set_size_inches(16,7)
    for j in range(len(MODS)):
        ax.plot(TIMESCONVRA_ALL[j],FULLCONV_RAINT_ALL[j])
    ax.set_ylabel('Ra internally')
    ax.set_xlabel('time (Myr)')
    ax.set_yscale('log')
    #plt.savefig('/Users/mkerr/VenusResearch/2022/sept2022_analysis/OCT4_l-Ra-slab_vs_time/Ra_eff.png',dpi=400)
    plt.show()
    plt.clf()
    plt.close()

    fig, ax = plt.subplots(1)
    fig.set_size_inches(16,7)
    for j in range(len(MODS)):
        ax.plot(TIMESCONVRA_ALL[j],FULLCONV_RACOMB_ALL[j])
    ax.set_ylabel('Ra combo (internal and bot)')
    ax.set_xlabel('time (Myr)')
    ax.set_yscale('log')
    #plt.savefig('/Users/mkerr/VenusResearch/2022/sept2022_analysis/OCT4_l-Ra-slab_vs_time/Ra_eff.png',dpi=400)
    plt.show()
    plt.clf()
    plt.close()

    ##############################################################################