import sys
import csv


STAGYY_WRITE_BL_RAC_PARS=int(sys.argv[2])

if STAGYY_WRITE_BL_RAC_PARS:
    STAGYY_OUTPUT_FOLDER='/Users/mkerr/VenusResearch/2023/STAGYY/mad_stag/_PUBMODS/TBL/'
    STAGYY_MODS_SPT1=[270,191,92,205,119,74,44] 
else:
    STAGYY_OUTPUT_FOLDER='/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/fromPerl/'
    STAGYY_MODS_SPT1=[234,234,234,234,234,234,234] #[843,1019,976,528,517,460,485] 
STAGYY_OUTPUT_MESH_FOLDER= '/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/'



STAGYY_OUTPUT_MODS=['1300_1/','1400_1/','1500_1/','1600_1/','1700_1/','1800_1/','1900_1/']
STAGYY_MODS_DIMS= ['512x128','512x128','512x128','1024x256','1024x256','1024x256','1024x256'] 
STAGYY_MODS_SPT0=[0,0,0,0,0,0,0] 
STAGYY_MODS_TCMB=[1400,1500,1600,1700,1800,1900,2000]
STAGYY_MODS_TM=[1300,1400,1500,1600,1700,1800,1900]

# FOUND FROM A HIGH-TIME-RES MODEL from t=0 to convection onset
CONV_ONSET_FILE='/Users/mkerr/VenusResearch/2023/MEASURE_KERR_2023/Deep-MPATK/output/conv_onset/'
STAGYY_WRITE_RA_EFF=int(sys.argv[3])
RA_EFF_FILE='/Users/mkerr/VenusResearch/2023/MEASURE_KERR_2023/Deep-MPATK/output/ra_eff/'
CONV_ONSET_MOD_PARS='/Users/mkerr/VenusResearch/2023/MEASURE_KERR_2023/Deep-MPATK/output/conv_onset/'

STAGYY_OUTPUT_MODS_0=0
STAGYY_OUTPUT_MODS_1=7


# TODO: index-ify read time and plot fields sphereical

# READ_TIME.PY 
READSTAGBIN_fields = ['eta','t'] # field name 't', 'eta', doesnt work with velocity/pressure yet.
READSTAGBIN_tstep0=0
READSTAGBIN_tstep1=40 # time-step
READSTAGBIN_tstep_step=1
READSTAGBIN_y_dim = 512 # grid dim on horz axis
READSTAGBIN_z_dim = 128 # grid dim on vert axis
READSTAGBIN_nz_cores = 2 #2 #number of core-blocks vertically
READSTAGBIN_ny_cores = 8 #8 # number of core-blocks horizontally
if STAGYY_WRITE_BL_RAC_PARS:
    READSTAGBIN_nz_cores = 1 #2 #number of core-blocks vertically
    READSTAGBIN_ny_cores = 1 #8 # number of core-blocks horizontally
READSTAGBIN_meshdir='/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/STAGmesh'+str(READSTAGBIN_y_dim)+'x'+str(READSTAGBIN_z_dim)+'.csv'
READSTAGBIN_showPlot_bool=0
READSTAGBIN_writetocsv_bool=1


# PLOT_FIELDS_SPHERICAL.PY
PFS_printVTR=0
PFS_savefigures=1
PFS_showfigures=0
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