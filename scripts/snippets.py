import set_DeepMPATK_pars as par
from matplotlib.cbook import flatten
from sklearn.cluster import KMeans
from sklearn import metrics as skmet;
import numpy as np
from statistics import mode
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import csv
import os
import stagyypythonmodule as spm


DIR=par.STAGYY_OUTPUT_FOLDER
TIMEDIR= par.STAGYY_OUTPUT_FOLDER

MESHDIR = par.STAGYY_OUTPUT_MESH_FOLDER

CSVWRITEDIR=par.SNIPPETS_CSVWRITEDIR
FIGSAVEDIR=par.SNIPPETS_FIGSAVEDIR

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
###############

################################
########  DEFINE MODULES ###### 
################################

MODS=par.STAGYY_OUTPUT_MODS
TCMB=par.STAGYY_MODS_TCMB
TM = par.STAGYY_MODS_TM
SPT0 = par.STAGYY_MODS_SPT0
SPT1 =  par.STAGYY_MODS_SPT1
BL_X_KM = par.STAGYY_MODS_BL_X_KM
BL_X_PT = par.STAGYY_MODS_BL_X_PT
RAC = par.STAGYY_MODS_RAC
DIMS = par.STAGYY_MODS_DIMS
TAU = par.STAGYY_MODS_TAU
N_ONSET=par.STAGYY_MODS_N_ONSET

if BL_X_KM==[]:
    raise ValueError('Cannot run snippets without getting the TBL, RA_C, TAU, and T_ONSET lists')

################################################################################################
################################################################################################

startindmod=par.STAGYY_OUTPUT_MODS_0
endin=par.STAGYY_OUTPUT_MODS_1
#startingTS=int(N_ONSET[startindmod])
#endingTS=SPT1[startindmod]
jump=par.SNIPPETS_JUMP

MODS=MODS[startindmod:endin]
TCMB=TCMB[startindmod:endin]
TM=TM[startindmod:endin]
SPT0=SPT0[startindmod:endin]
SPT1=SPT1[startindmod:endin]
RAC=RAC[startindmod:endin]
BL_X_KM=BL_X_KM[startindmod:endin]
BL_X_PT=BL_X_PT[startindmod:endin]
DIMS=DIMS[startindmod:endin]
TAU=TAU[startindmod:endin]

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
# LOOP THROUGH DIRECTORIES
##############################################################################
seg=0

for i in range(len(DIRS)): 
    file = open(CSVWRITEDIR+MODS[i][:-1]+'.txt', 'w')
    file.close()

HD_OV_CON_SZ=[]
HD_OV_CON_NO=[]
TM_SNP=[]
AVG_PL_NO=[]

# ******** ------- ********** ----------- ************* -----------------
# ******** ------- ********** ----------- ************* -----------------
for i in range(len(DIRS)): #for each of the models
    dir=DIRS[i]
    print('******** CURRENT DIR:',dir,' *********')
    BLTHICKkm_itpl=[]
    BLTHICKp_itpl=[]
    TIMES=[]
    RA=[]
    RA_smooth=[]
    RA_BIN=[]
    PLMNO=[]
    PLMTK=[]
    NO_TK_plumes=[]
    seg=0
    startingTS=int(N_ONSET[i])
    endingTS=SPT1[i]
    for snip_ind in range(startingTS,endingTS,jump):
        print('SEG:', seg)
        spt0=snip_ind
        spt1=snip_ind+jump

        if spt1>endingTS:
            spt1=endingTS

        print('snippet:',snip_ind,': ', spt0,'->',spt1)

        file = open(CSVWRITEDIR+MODS[i][:-1]+'.txt', 'a')
        file.write('\n snippet:'+str(seg)+': '+ str(spt0)+'->'+str(spt1)+'\n')
        nums=NUMS[i][spt0:spt1]
        print('len nums:',len(nums))
        t0 = TM[i]
        tcmb=TCMB[i]

        Times=spm.getTimeARR(TDIRS[i],spt0,spt1,1)  #time directory, step size
        print('SPT1:', spt1)
        len_Time=len(Times)
        print(len_Time)

        Ra_c=RAC[i]
        blX_km=BL_X_KM[i]
        blX_pt=BL_X_PT[i]
        print(blX_km,blX_pt)
        meshdim=DIMS[i]
        if meshdim=='1024x256':
            mesh_y_dim=1024
        elif meshdim=='512x128':
            mesh_y_dim=512
        elif meshdim=='256x64':
            mesh_y_dim=256
        else:
            print('MESH ERROR')
            break


        #OPEN MESH FILE	
        Tablemesh=[]
        with open('/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/STAGmesh'+meshdim+'.csv') as f:
            reader = csv.reader(f)
            for col in reader:
                Tablemesh.append(col)
                        
        print(len(Tablemesh))
        print(len(Tablemesh[0]))
        y=[float(i) for i in Tablemesh[0][0:mesh_y_dim]] # get the radii coords
        z=[float(i) for i in Tablemesh[0][mesh_y_dim::2*mesh_y_dim]] # get the z coords
        Y=len(y)
        Z=len(z)
        arrY=[i for i in range(1,Y+1)]
        arrZ=[i for i in range(1,Z+1)]
        eighth_of_dom=(mesh_y_dim//4)//8
        mm0=2*eighth_of_dom
        mm1=5*eighth_of_dom
        
        KMAXNUM=par.SNIPPETS_KMAXNUM
        WALLSIZE= int(np.floor(1.5*blX_pt)) #mesh_y_dim//128 
        WALLSIZE= int(np.floor(2*blX_pt)) # FOR 1300 ONLY
        WALLBOUNDL=WALLSIZE
        WALLBOUNDR=-1*WALLSIZE
        #############################################################
        midmantlen = [mm0,mm1]
        dims = [Y,Z]
        # (DEFINED ABOVE) geom_pars = [domain_THICK,air_THICK,core_THICK,D_MantlenStag,mantle_THICK]
        # (DEFINED ABOVE) phys_pars=[a_exp,g,Cp,kappa,rho,therm_conductivity,H_flux,H_mass,ra_num_const]
        pars = [dims,phys_pars,geom_pars,midmantlen]
        RAD = spm.getRAD(TIMEDIR+MODS[0],pars)

        BLTHICKkm_1mod_itpl=[]
        BLTHICKp_1mod_itpl=[]

        TIMES.append(Times)
        RAplot=[]
        RAcritBIN=[]
        RAgrpPLT=[]

        PLMNO_1mod=[]
        PLMTK_1mod=[]
        NO_TK_plumes_1mod=[]

        BLdim_p=0
        BLdim_km=0.1

        #################################
        for n in nums: #for each time-step of each model
            N=int(n)
            #print('TIMESTEP:  ',n) # PRINT STATEMENT
            #writerAVG.writerow([t]) # add the model time to the output csv file

            #################################
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
            TableV=[];	
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
            delt=tcmb-TMEAN

            ########################################

            THIGH=TMEAN+(0.1*delt)                                  # % TBL TEMP VALUE
            ### FROM TURCOTTE AND SHUBERT PG 184 ON HEAT TRANSFER
         # ******************************************************************************************************************************************************

            BLDIM_km=blX_km
            BLDIM_p=blX_pt
         # ******************************************************************************************************************************************************


            RAloc_2 = spm.getRAlslice(BLDIM_p,BLDIM_km,delt,TableV,pars) #polynomial
            RAloc_2=RAloc_2[WALLBOUNDL:WALLBOUNDR]
            RAplot.append(RAloc_2)
            Ra_crit=Ra_c
            RA_binSLC=[el>=Ra_crit for el in RAloc_2]
            RAcritBIN.append(RA_binSLC)  # change here for plume vs slab detection

         # ******************************************************************************************************************************************************

            #[pnum,pwid]=plumeNo_Th(Ra_crit,BLkm_itpl,RAloc_2,1)
            [pnum,pwid]=spm.plumeNo_Th(Ra_crit,BLDIM_km,RAloc_2,1,pars)
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

            ##############################################################################
            # PLOTS (each Time step)
            ##############################################################################
            if par.SNIPPETS_plotRalocandRac: # plot the boundary layer thickness across the Y domain
                fig, ax = plt.subplots(1)
                fig.set_size_inches(5,8)
                fig.suptitle(MODS[i][:-1])
                ax.plot(WALLBOUNDL+np.arange(len(RAloc_2)),RAloc_2,color='k')
                ax.plot(np.arange(Y), np.ones(Y)*Ra_c,color='red')
                plt.show()
                plt.clf()
                plt.close()

            ##############################################################################

        RA.append(RAplot)
        RA_BIN.append(RAcritBIN)
        ####################################
        print('TBL THICKNESS:',np.mean(BLDIM_km)) # PRINT STATEMENT
        BLTHICKkm_itpl.append(BLTHICKkm_1mod_itpl)
        BLTHICKp_itpl.append(BLTHICKp_1mod_itpl)

        NO_TK_plumes.append(NO_TK_plumes_1mod)
        PLMNO.append(PLMNO_1mod)
        allplumes_1mod=sorted(flatten(PLMTK_1mod))
        allplumes=np.array(allplumes_1mod)
        allplumes.reshape(-1,1)
        
        PLMTK.append(sorted(allplumes))


        #  * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        # * * * * * * * *     KMEANS      * * * * * * * * * * * * 
        # * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

        kmax=KMAXNUM
        CENTROIDS_1mod=[]
        print('NUMBER OF PLUMES:',len(allplumes))
        data = allplumes.reshape(-1,1)
        stats_empty=1

        if len(allplumes)<kmax:
            continue

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
                if origDisp==0.0:
                    gap = np.log(np.mean(refDisps))
                else:
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
                file.write('['+str(seg)+','+str(ii)+','+str([x[0] for x in clumpdata])+','+str([x[1] for x in clumpdata])+','+str([x[2] for x in clumpdata])+','+str([Times[0], Times[-1]])+']\n')
                if ii==2 and par.SNIPPETS_makeK2stats:
                    sizes=[x[0] for x in clumpdata]
                    numbs=[x[2] for x in clumpdata]
                    HD_OV_CON_SZ.append(sizes[1]/sizes[0])
                    HD_OV_CON_NO.append(numbs[1]/numbs[0])
                    TM_SNP.append(Times[0]/TAU[i])
                    AVG_PL_NO.append(np.sum(numbs)/len_Time)
                    print(TM_SNP)
                    print(HD_OV_CON_SZ)
                    print(HD_OV_CON_NO)
                    print(AVG_PL_NO)
                if ii>1 and par.SNIPPETS_showClustergroupsforK:
                    fig, ax = plt.subplots(1)
                    fig.set_size_inches(5,3)
                    fig.suptitle(MODS[i][:-1])
                    ax.scatter([x[0] for x in clumpdata],[x[2] for x in clumpdata])
                    ax.errorbar([x[0] for x in clumpdata],[x[2] for x in clumpdata],xerr=[x[1] for x in clumpdata])
                    plt.show()
                    plt.clf()
                    plt.close()
                if origDisp==0.0:
                    break

            SSdiff_elbow=np.array(SSdiff_elbow)
        else:
            print('No plumes')
            SSdiff_CC=np.arange(1,kmax,1)


         # ******************************************************************************************************************************************************

        arrY1=arrY
        arrY=arrY[WALLBOUNDL:WALLBOUNDR+1]
        if len(CENTROIDS_1mod)!=0:
            testR=1
            for cc in range(1,len(CENTROIDS_1mod)+1,1):
                RAgrpPLT=[]
                for tt in range(len(RAcritBIN)): #JAN5
                    intMask=[int(x) for x in RAcritBIN[tt]]
                    newbinSLC,Y_HEAD_COORDS=spm.clumpedSlice(intMask,NO_TK_plumes_1mod[tt][0],NO_TK_plumes_1mod[tt][1],CENTROIDS_1mod[cc-1][0],CENTROIDS_1mod[cc-1][1])
                    RAgrpPLT.append(newbinSLC)
                RApts=RAplot #[spt0:spt1+1] # SPLIT
                fileName=MODS[i][:-1]
                fig, ax = plt.subplots(1) # JAN4
                fig.set_size_inches(16,6)
                fig.suptitle(fileName+" k="+str(cc)+'   seg='+str(seg))
                #print(len(arrY),len(TIMES[i][spt0:spt1]),len(RApts),len(RApts[0]))
                XCMB,TIMEevo = np.meshgrid(arrY,Times)
                ax.pcolor(XCMB,TIMEevo,RApts,cmap='Blues',vmin=0,vmax=5000) #SPLIT
                c=ax.pcolor(XCMB,TIMEevo,RAgrpPLT,cmap=cmap_rb,vmin=0,vmax=cc) # groups #SPLIT
                #ax.set_ylim([TIMES[i][spt0]-1,TIMES[i][spt1]+1])
                cax = plt.axes([0.93, 0.1, 0.02, 0.8])
                plt.savefig(FIGSAVEDIR+MODS[i][:-1]+'/seg'+str(seg)+'_k'+str(cc)+'.png',dpi=300)
                #plt.show()
                plt.clf()
                plt.close()
        arrY=arrY1
        # ******************************************************************************************************************************************************

        ##############################################################################
        # PLOTS (each model)
        ##############################################################################
        plotss=1
        if plotss==1:
            fig, ax = plt.subplots(6)
            fig.set_size_inches(5,7)
            fig.suptitle(MODS[i][:-1])
            numkplot=len(SSdiff_CC)+1
            ax[0].plot(np.arange(1,numkplot,1),SSdiff_CC)
            ax[1].plot(np.arange(1,numkplot,1),SSdiff_CC)
            ax[2].plot(np.arange(1,numkplot,1),DB_1mod)
            ax[3].plot(np.arange(1,numkplot,1),CH_1mod)
            ax[4].plot(np.arange(1,numkplot,1),silo_1mod)
            ax[5].plot(np.arange(1,numkplot,1),GAP_1mod)
            ax[1].set_yscale('log')
            ax[0].set_ylabel('ss diff')
            ax[1].set_ylabel('ss diff -log')
            ax[2].set_ylabel('DB')
            ax[3].set_ylabel('CH')
            ax[4].set_ylabel('silohuette')
            ax[5].set_ylabel('gapstat')
            ax[0].set_xlabel('k')
            plt.savefig(FIGSAVEDIR+MODS[i][:-1]+'/kOPT_seg'+str(seg)+'_k'+str(cc)+'.png')
            #plt.show()
            plt.clf()
            plt.close()
        file.close()
        seg=seg+1
    # ******** ------- ********** ----------- ************* -----------------
# ******** ------- ********** ----------- ************* -----------------

fig, ax = plt.subplots(2)
fig.set_size_inches(4,7)
ax[0].plot(HD_OV_CON_NO)
ax[1].plot(HD_OV_CON_SZ)

plt.show()
plt.clf()
plt.close()