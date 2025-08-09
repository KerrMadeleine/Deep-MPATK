import sys
import csv

"""
This parameter file takes in the flag inputs to the python scripts to determine which "route" down the 
pipeline the code base will take. STAGYY_WRITE_BL_RAC_PARS is better named "use HiRes models" since
it is the high-resolution-in-time models that are interpolated in time for the critical boundary layer
thickness and the critical Rayleigh number of that model.

This first global variable "STAGYY_WRITE_BL_RAC_PARS" is a switch (ON=use high-res-in-time models, 
OFF = use regular, full time models).

SPT is short for snippet and SPT0 is the output file number for the beginning of the snippet
and SPT1 is the output file number for the end of the snippet. A snippet is an interval of timesteps
between the initial and final times. 

MODS are models, DIMS are dimensions, TCMB is the temperature at the core-mantle boundary condition
TM is the initialized mantle potential temperature.
TODO: more description of global variable names
"""


STAGYY_WRITE_BL_RAC_PARS=int(sys.argv[2])

if STAGYY_WRITE_BL_RAC_PARS:
    #STAGYY_OUTPUT_FOLDER='/Users/mkerr/VenusResearch/2023/STAGYY/mad_stag/_PUBMODS/TBL/'
    STAGYY_OUTPUT_FOLDER='/Volumes/easystore/ProjectFolders/STAGYY/mad_stag/_PUBMODS/TBL/'
    STAGYY_MODS_SPT1=[270,191,92,205,119,74,44] 
else:
    #STAGYY_OUTPUT_FOLDER='/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/fromPerl/'
    STAGYY_OUTPUT_FOLDER='/Volumes/easystore/ProjectFolders/STAGYY/Analysis/fromPerl/'
    STAGYY_MODS_SPT1=[843,1019,976,528,517,460,485] 

#STAGYY_OUTPUT_MESH_FOLDER= '/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/'
STAGYY_OUTPUT_MESH_FOLDER= '/Volumes/easystore/ProjectFolders/STAGYY/Analysis/'

"""
The global variables defined here 
"""

STAGYY_OUTPUT_MODS=['1300_1/','1400_1/','1500_1/','1600_1/','1700_1/','1800_1/','1900_1/']
STAGYY_MODS_DIMS= ['512x128','512x128','512x128','1024x256','1024x256','1024x256','1024x256'] 
STAGYY_MODS_SPT0=[842,1019,1922,601,685,631,666] 
#STAGYY_MODS_SPT0=[842,1019,976,528,517,460,485] 
STAGYY_MODS_TCMB=[1400,1500,1600,1700,1800,1900,2000]
STAGYY_MODS_TM=[1300,1400,1500,1600,1700,1800,1900]
STAGYY_DEEPMPATK_OUTPUT_DIR='/Volumes/easystore/ProjectFolders/Deep-MPATK/output/'

DOMAIN_THICKNESS = 2942.e3 # m
AIR_THICKNESS = 0 # thickness of sticky-air if needed
CORE_RADIUS = 3110.e3 # m 

THERM_EXPANS = 3.e-5 # K^-1 , alpha
GRAV_ACCEL=8.87 # m/s^2 , g
SPECF_HEAT_CONST_PRESS=1200 # J/K/kg , Cp
THERM_DIFFUS=1.e-6 #m^2/s , kappa
REF_DENSITY=3300 # kg/m^3 , rho
THERM_COND = THERM_DIFFUS*REF_DENSITY*SPECF_HEAT_CONST_PRESS
HEAT_PROD_FLUX = 0 #1.623e-8 # W/m^3
HEAT_PROD_MASS = 0 #4.9185e-12 # W/kg

MANTLE_THICKNESS=DOMAIN_THICKNESS-AIR_THICKNESS
ra_num_const = (MANTLE_THICKNESS)**3 * REF_DENSITY**2 * SPECF_HEAT_CONST_PRESS*GRAV_ACCEL*THERM_EXPANS/THERM_COND
MANTLE_THICKNESS_KM=MANTLE_THICKNESS/1000

# FOUND FROM A HIGH-TIME-RES MODEL from t=0 to convection onset
#CONV_ONSET_FILE='/Users/mkerr/VenusResearch/2023/MEASURE_KERR_2023/Deep-MPATK/output/conv_onset/'
CONV_ONSET_FILE='/Volumes/easystore/ProjectFolders/Deep-MPATK/output/conv_onset/'
STAGYY_WRITE_RA_EFF=int(sys.argv[3])
#RA_EFF_FILE='/Users/mkerr/VenusResearch/2023/MEASURE_KERR_2023/Deep-MPATK/output/ra_eff/'
RA_EFF_FILE='/Volumes/easystore/ProjectFolders/Deep-MPATK/output/ra_eff/'
#CONV_ONSET_MOD_PARS='/Users/mkerr/VenusResearch/2023/MEASURE_KERR_2023/Deep-MPATK/output/conv_onset/'
CONV_ONSET_MOD_PARS='/Volumes/easystore/ProjectFolders/Deep-MPATK/output/conv_onset/'



STAGYY_OUTPUT_MODS_0=3
STAGYY_OUTPUT_MODS_1=4

# TODO: index-ify read time and plot fields sphereical

# READ_TIME.PY 
READSTAGBIN_fields = ['eta','t'] # field name 't', 'eta', doesnt work with velocity/pressure yet.
READSTAGBIN_tstep0= 600 #200
READSTAGBIN_tstep1= 602# time-step
READSTAGBIN_tstep_step=1
READSTAGBIN_y_dim = 1024 # grid dim on horz axis
READSTAGBIN_z_dim = 256 # grid dim on vert axis


"""
For reading the StagYY binary files, it is important to know how the domain is split between cores.
1 x 1 core geometry for all hi-time-res models  
2 x 8 for 1300 - 1500 K all times x
2 x 16 for 1600 K all times x (602)
2 x 8 for 1700 K, initially (0-356), then 2 x 16 (357 to 512) x | o  (685)
2 x 8 for 1800 K, initially (0-283), then 2 x 16 ( 283 to 460) x | o (631)
2 x 8 for 1900 K all times o (666)

Important to know that some models were run with 16 cores, then 32 so everythings a bit wonky here
"""
READSTAGBIN_nz_cores = 2 #2 #2 #number of core-blocks vertically
READSTAGBIN_ny_cores = 16 #8 #8 # number of core-blocks horizontally
if STAGYY_WRITE_BL_RAC_PARS:
    # The hi res models were run locally on my laptop, so the structure of the binary files 
    # is different based on the way the domain is divided for each CPU
    READSTAGBIN_nz_cores = 1 #number of core-blocks vertically
    READSTAGBIN_ny_cores = 1 # number of core-blocks horizontally
READSTAGBIN_meshdir='/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/STAGmesh'+str(READSTAGBIN_y_dim)+'x'+str(READSTAGBIN_z_dim)+'.csv'
READSTAGBIN_showPlot_bool=0
READSTAGBIN_writetocsv_bool=1




# PLOT_FIELDS_SPHERICAL.PY
PFS_printVTR=0
PFS_savefigures=0
PFS_showfigures=1
PFS_diffFromMeanT=1





# READ_RPROF.PY
READ_RPROF_SHOWFIG=0
READ_RPROF_SAVEFIG=0





# PIPELINE.PY
PIPELINE_HAS_BLX_RA=int(sys.argv[1])
PIPELINE_WRITE_CONV_REG=STAGYY_WRITE_RA_EFF
print(STAGYY_WRITE_RA_EFF,PIPELINE_WRITE_CONV_REG)

if PIPELINE_HAS_BLX_RA:
    PIPELINE_calcKmeansBool=1
    PIPELINE_CALC_CONV_REG=1
    PIPELINE_startaftertonset=1
    PIPELINE_makeBLplots=0
    PIPELINE_computecoldTBL=0
    PIPELINE_wallsizefactor=1
    PIPELINE_plotRaloc_w_mask=0
    PIPELINE_tempKellogg=0
    PIPELINE_plotKmeansKellogg=0

    PAR_DATA=[]
    with open(CONV_ONSET_FILE+'P1_pars.txt') as f:
        reader = csv.reader(f,delimiter='\n')
        for col in reader:
            PAR_DATA.append(eval(col[0]))

    STAGYY_MODS_BL_X_KM= PAR_DATA[0]   #[346.4,203.1]          # [1746.75,1026.12,616.55,346.4,203.1,125.063,83.3255,60.472,45.555]
    STAGYY_MODS_BL_X_PT= PAR_DATA[1]   #[19.0016,11.079]       # [43.63,27.67,33.392,19.0016,11.079,13.855,9.0989,6.48,4.765]
    STAGYY_MODS_RAC= PAR_DATA[2]       #[1703,1888]            # [1107.01,1253.2,1397.77,1703,1888,1940.05,2133.89,2627.1,3197.71]
    STAGYY_MODS_N_ONSET= PAR_DATA[3]    #[10,8]                 # [0,223,37,10,8,7,11,50,45]
    STAGYY_MODS_TAU= PAR_DATA[4]      #[730,230]              # [41460e3,2.773e3,2400,730,230,82,37,15,8.6]
    # print statments to sanity-check/debug
    print(STAGYY_MODS_BL_X_KM)
    print(STAGYY_MODS_BL_X_PT)
    print(STAGYY_MODS_RAC)
    print(STAGYY_MODS_N_ONSET)
    print(STAGYY_MODS_TAU)

else:
    PIPELINE_calcKmeansBool=0
    PIPELINE_CALC_CONV_REG=1
    PIPELINE_makeBLplots=1
    PIPELINE_computecoldTBL=1
    PIPELINE_startaftertonset=0
    PIPELINE_wallsizefactor=0
    PIPELINE_plotRaloc_w_mask=0
    PIPELINE_tempKellogg=0
    PIPELINE_plotKmeansKellogg=0

# parameters
PIPELINE_KMAXNUM=10
PIPELINE_BLplotfreq = 2000
PIPELINE_TprofFreq = 2000
PIPELINE_VprofFreq = 2000
PIPELINE_tempKellogg_minusMean=1

# misc plotting
PIPELINE_showStatsandScores=0
PIPLELINE_makeRavTimePlot=0
PIPELINE_makeBLvsTimeplot=0
PIPELINE_makeannotatedTprof=0
PIPELINE_showViscprof=0
PIPELINE_plotallmodsRavTime=0




# SNIPPETS.PY
SNIPPETS_CSVWRITEDIR='/Users/mkerr/VenusResearch/2023/MEASURE_KERR_2023/Deep-MPATK/output/snippets/csvFiles/'
SNIPPETS_FIGSAVEDIR='/Users/mkerr/VenusResearch/2023/MEASURE_KERR_2023/Deep-MPATK/output/snippets/figures/'

SNIPPETS_JUMP = 50
SNIPPETS_HAS_BL_RA_TAU=0

SNIPPETS_makeK2stats = 0
SNIPPETS_plotRalocandRac=0
SNIPPETS_showClustergroupsforK=0
SNIPPETS_KMAXNUM = 6