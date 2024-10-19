import set_DeepMPATK_pars as par
from sklearn.cluster import KMeans
from sklearn import metrics as skmet;
import numpy as np
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import csv
import os
import stagyypythonmodule as spm
from cmcrameri import cm

"""
This script plots the temperature, viscosity and density fields over the domain
and saves these figures.
"""

DIR=par.STAGYY_OUTPUT_FOLDER
TIMEDIR= par.STAGYY_OUTPUT_FOLDER
MESHDIR = par.STAGYY_OUTPUT_MESH_FOLDER

# TODO: move this to global parameters
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

DIR=par.STAGYY_OUTPUT_FOLDER
TIMEDIR= par.STAGYY_OUTPUT_FOLDER
MODS=par.STAGYY_OUTPUT_MODS
SPT0 = par.STAGYY_MODS_SPT0
SPT1 = par.STAGYY_MODS_SPT1
DIMS = par.STAGYY_MODS_DIMS

savefigures=par.PFS_savefigures
showfigures=par.PFS_showfigures
printVTR=par.PFS_printVTR
startindmod=par.STAGYY_OUTPUT_MODS_0
endin=par.STAGYY_OUTPUT_MODS_1
fields=par.READSTAGBIN_fields
yesRHO=1
if len(fields)==2:
    yesRHO=0



####################################################################################
MODS=MODS[startindmod:endin]
SPT0=SPT0[startindmod:endin]
SPT1=SPT1[startindmod:endin]
DIMS=DIMS[startindmod:endin]

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
os.makedirs(par.STAGYY_DEEPMPATK_OUTPUT_DIR+'figures', exist_ok = True)
for i in range(len(DIRS)): #for each of the models
    print('Making plots for:', MODS[i])
    save_fig_dir=par.STAGYY_DEEPMPATK_OUTPUT_DIR+'figures/'+MODS[i][:-1]
    os.makedirs(save_fig_dir, exist_ok = True)

    dir=DIRS[i]
    nums=NUMS[i]
    spt0=SPT0[i]
    spt1=SPT1[i]

    Times=spm.getTimeARR(TDIRS[i],spt0,spt1,1)  #time directory, step size
    nums_int = [int(n)for n in nums[spt0:]]
    meshdim=DIMS[i]
    print("Mesh dimensions:", meshdim)
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
    with open(MESHDIR+'STAGmesh'+meshdim+'.csv') as f:
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


    ##############################################################################
    midmantlen = [mm0,mm1]
    dims = [Y,Z]
    # (DEFINED ABOVE) geom_pars = [domain_THICK,air_THICK,core_THICK,D_MantlenStag,mantle_THICK]
    # (DEFINED ABOVE) phys_pars=[a_exp,g,Cp,kappa,rho,therm_conductivity,H_flux,H_mass,ra_num_const]
    pars = [dims,phys_pars,geom_pars,midmantlen]
    RAD = spm.getRAD(TIMEDIR+MODS[0],pars)

    print('******** CURRENT DIR:',dir,' *********') # PRINT STATEMENT
    arrY=np.array(arrY)*(np.pi/mesh_y_dim)

    ###########################################################################
    for n in nums[spt0:spt1+1]: #for each time-step of each model
        print('TIMESTEP:  ',n,end='\r') # PRINT STATEMENT
        N=int(n)
        #writerAVG.writerow([t]) # add the model time to the output csv file

        #################################
        # READ THE MOD/TS FILE
        openFileT= dir+'t_'+n+'.csv' #temp file to open
        openFileV= dir+'eta_'+n+'.csv' #visc file to open

        if yesRHO:
            openFileR= dir+'rho_'+n+'.csv' #visc file to open

        #OPEN TEMP FILE
        TableT=[]
        with open(openFileT) as f:
            reader = csv.reader(f,delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)
            for col in reader:
                #print(col)
                TableT.append(col)
                
        #OPEN VISC FILE	
        TableV=[];	
        with open(openFileV) as f:
            reader = csv.reader(f, delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)
            for col in reader:
                #print(col)
                TableV.append(col)

        if yesRHO:
            TableR=[];	
            with open(openFileR) as f:
                reader = csv.reader(f, delimiter=' ',quoting=csv.QUOTE_NONNUMERIC)
                for col in reader:
                    #print(col)
                    TableR.append(col)

        #print("shape of table T",len(TableT),len(TableT[0]))
        TableT=np.ravel(TableT,order='C')
        TableV=np.ravel(TableV,order='C')

        if yesRHO:
            TableR=np.ravel(TableR,order='C')

        TGRAD=[spm.trimmedMean(spm.getRadSlice(i,TableT,pars)) for i in range(Z)]
        if yesRHO:
            RGRAD=[spm.trimmedMean(spm.getRadSlice(i,TableR,pars)) for i in range(Z)]
        VGRAD=[spm.trimmedMeanV(spm.getRadSlice(i,TableV,pars)) for i in range(Z)]
        if printVTR:
            print(VGRAD)
            print(TGRAD)
            if yesRHO:
                print(RGRAD)

        TMEAN = np.mean(TGRAD[mm0:mm1])

        RAD_PLT=np.array(RAD)+3110

        if yesRHO:
            R=np.array([spm.getRadSlice(i,TableR,pars) for i in range(Z)])
            R=np.fliplr(R)
        
        T=np.array([spm.getRadSlice(i,TableT,pars) for i in range(Z)])
        V=np.log10(np.array([spm.getRadSlice(i,TableV,pars) for i in range(Z)])) #-np.log10(VMEAN)
        V=np.fliplr(V)
        T=np.fliplr(T)
        

        fig, ax = plt.subplots(1,subplot_kw={'projection': 'polar'})
        fig.set_size_inches(9.5,7)
        fig.suptitle('T at t='+str(Times[N-spt0])[:6]+' Myr',size=20)
        XCMB,YDEP = np.meshgrid(arrY,RAD_PLT)
        if par.PFS_diffFromMeanT:
            c=ax.pcolor(XCMB,YDEP,T-TMEAN,cmap=cm.lajolla,vmin=-200,vmax=200)
        else:
            c=ax.pcolor(XCMB,YDEP,T,cmap=cm.lajolla,vmin=1000,vmax=3500)
        ax.set_rticks([3110,3110+mantle_THICK])  # Less radial ticks
        ax.set_thetagrids([])
        ax.set_thetalim(0,np.pi)
        cax = plt.axes([0.9, 0.25, 0.02, 0.5])
        fig.colorbar(c,cax=cax)
        plt.subplots_adjust(left=0.07, bottom=0, right=0.9, top=1, wspace=0, hspace=0)
        if savefigures:
            plt.savefig(save_fig_dir+'/'+str(N)+'_T.png',dpi=400)
        if showfigures:
            plt.show()
        plt.clf()
        plt.close()

        fig, ax = plt.subplots(1,subplot_kw={'projection': 'polar'})
        fig.set_size_inches(9.5,7)
        fig.suptitle('V at t='+str(Times[N-spt0])[:6]+' Myr',size=20)
        XCMB,YDEP = np.meshgrid(arrY,RAD_PLT)
        c=ax.pcolor(XCMB,YDEP,V,cmap=cm.cork_r,vmin=18,vmax=25)
        ax.set_rticks([3110,3110+mantle_THICK])  # Less radial ticks
        ax.set_thetagrids([])
        ax.set_thetalim(0,np.pi)
        cax = plt.axes([0.9, 0.25, 0.02, 0.5])
        fig.colorbar(c,cax=cax)
        plt.subplots_adjust(left=0.07, bottom=0, right=0.9, top=1, wspace=0, hspace=0)
        if savefigures:
            plt.savefig(save_fig_dir+'/'+str(N)+'_V.png',dpi=400)
        if showfigures:
            plt.show()
        plt.clf()
        plt.close()

        if yesRHO:
            fig, ax = plt.subplots(1,subplot_kw={'projection': 'polar'})
            fig.set_size_inches(9.5,7)
            fig.suptitle('rho at t='+str(Times[N-spt0])[:6]+' Myr',size=20)
            XCMB,YDEP = np.meshgrid(arrY,RAD_PLT)
            c=ax.pcolor(XCMB,YDEP,R,cmap=cm.cork_r,vmin=3000,vmax=3400)
            ax.set_rticks([3110,3110+mantle_THICK])  # Less radial ticks
            ax.set_thetagrids([])
            ax.set_thetalim(0,np.pi)
            cax = plt.axes([0.9, 0.25, 0.02, 0.5])
            fig.colorbar(c,cax=cax)
            plt.subplots_adjust(left=0.07, bottom=0, right=0.9, top=1, wspace=0, hspace=0)
            if savefigures:
                plt.savefig(save_fig_dir+'/'+'_R.png',dpi=400)
            if showfigures:
                plt.show()
            plt.clf()
            plt.close()