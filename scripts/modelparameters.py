
# Global variables
DOMAIN_THICKNESS = 2942.e3 # m
TSURF = 740
AIR_THICKNESS = 0 # thickness of sticky-air if needed
CORE_RADIUS = 3110.e3 # m 
THERM_EXPANS = 3.e-5 # K^-1 , alpha
GRAV_ACCEL=8.87 # m/s^2 , g
SPECF_HEAT_CONST_PRESS=1200 # J/K/kg , Cp
THERM_DIFFUS=1.e-6 #m^2/s , kappa
REF_DENSITY=3300 # kg/m^3 , rho
THERM_COND = THERM_DIFFUS*REF_DENSITY*SPECF_HEAT_CONST_PRESS

MESH_DIR='/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/'
ANALYSIS_OUTPUT_DIR='/Volumes/easystore/ProjectFolders/Deep-MPATK/output/'
STAGYY_OUTPUT_FOLDER='/Volumes/easystore/ProjectFolders/STAGYY/Analysis/fromPerl/'
STAGYY_HIRES_OUTPUT_FOLDER='/Volumes/easystore/ProjectFolders/STAGYY/mad_stag/_PUBMODS/TBL/'
# CONV_ONSET_FILE='/Volumes/easystore/ProjectFolders/Deep-MPATK/output/conv_onset/'
# RA_EFF_FILE='/Volumes/easystore/ProjectFolders/Deep-MPATK/output/ra_eff/'
# CONV_ONSET_MOD_PARS='/Volumes/easystore/ProjectFolders/Deep-MPATK/output/conv_onset/'


# Model parameters
STAGYY_MODS=['1300_1/','1400_1/','1500_1/','1600_1/','1700_1/','1800_1/','1900_1/']
STAGYY_MODS_DIMS_STR = ['512x128','512x128','512x128','1024x256','1024x256','1024x256','1024x256'] 
STAGYY_MODS_YDIM = [512,512,512,1024,1024,1024,1024]
STAGYY_MODS_ZDIM = [128,128,128,256,256,256,256]
STAGYY_MODS_NMAX=[842,1019,1922,601,685,631,666] 
STAGYY_MODS_TCMB=[1400,1500,1600,1700,1800,1900,2000]
STAGYY_MODS_TM=[1300,1400,1500,1600,1700,1800,1900]

STAGYY_MODS_0_NZ = [2,2,2,2,2,2,2]
STAGYY_MODS_0_NY = [8,8,8,16,8,8,8]
STAGYY_MODS_NCORE_CHANGE=[0,0,0,0,356,283,0] 
STAGYY_MODS_1_NZ = [2,2,2,2,2,2,2]
STAGYY_MODS_1_NY = [8,8,8,16,16,16,8]

# For reading the StagYY binary files, it is important to know how the domain is split between cores.
# 1 x 1 core geometry for all hi-time-res models
# 2 x 8 for 1300 - 1500 K all times
# 2 x 16 for 1600 K all times
# 2 x 8 for 1700 K, initially (0-355), then 2 x 16 (356 to 512)
# 2 x 8 for 1800 K, initially (0-282), then 2 x 16 ( 283 to 460)
# 2 x 8 for 1900 K all times
# nz x ny

# Hires model parameters
STAGYY_MODS_HIRES_DIR='/Volumes/easystore/ProjectFolders/STAGYY/mad_stag/_PUBMODS/TBL/'
STAGYY_MODS_HIRES=['1300_1/','1400_1/','1500_1/','1600_1/','1700_1/','1800_1/','1900_1/']
STAGYY_MODS_HIRES_TCMB=[1400,1500,1600,1700,1800,1900,2000]
STAGYY_MODS_HIRES_TM=[1300,1400,1500,1600,1700,1800,1900]
STAGYY_MODS_HIRES_NMAX=[270,191,92,205,119,74,44]
STAGYY_MODS_HIRES_DIMS = ['512x128','512x128','512x128','1024x256','1024x256','1024x256','1024x256']
STAGYY_MODS_HIRES_YDIM = [512,512,512,1024,1024,1024,1024]
STAGYY_MODS_HIRES_ZDIM = [128,128,128,256,256,256,256]
STAGYY_MODS_CORE_NZ = [1,1,1,1,1,1,1]
STAGYY_MODS_CIRE_NY = [1,1,1,1,1,1,1]