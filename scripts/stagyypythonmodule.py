import numpy as np
import matplotlib.pyplot as plt
import csv


'''
This document has the functions used throughout the pipeline.
TODO: identify more functions that are built-in in numpy and replace the version i made myself
'''
#midmantlen = [mm0,mm1]
#dims = [Y,Z]
# (DEFINED ABOVE) geom_pars = [domain_THICK,air_THICK,core_THICK,D_MantlenStag,mantle_THICK]
# (DEFINED ABOVE) phys_pars=[a_exp,g,Cp,kappa,rho,therm_conductivity,H_flux,H_mass,ra_num_const]
#pars = [dims,phys_pars,geom_pars,midmantlen]


##############################################################################
# GLOBAL FUNCTIONS
##############################################################################
def trimmedMean(list):

	mean = np.mean(list)
	sd = np.std(list)
	final_list = [x for x in list if (x >= mean - sd)]
	final_list = [x for x in final_list if (x <= mean + sd)]

	mean=np.mean(final_list)
	sd = np.std(final_list)
	final_list2=[x for x in final_list if (x >= mean - sd)]
	final_list2 = [x for x in final_list2 if (x <= mean + sd)]

	mean=np.mean(final_list2)
	sd = np.std(final_list2)
	final_list3=[x for x in final_list2 if (x >= mean - sd)]

	
	return np.mean(final_list3)
	
def trimmedMeanV(listV):

    loglist=[np.log(i) for i in listV]
    mean = np.mean(loglist)
    sd = np.std(loglist)
    final_list = [x for x in loglist if (x >= mean - sd)]
    final_list = [x for x in final_list if (x <= mean + sd)]

    mean=np.mean(final_list)
    sd = np.std(final_list)
    final_list2=[x for x in final_list if (x >= mean - sd)]
    final_list2 = [x for x in final_list2 if (x <= mean + sd)]

    mean=np.mean(final_list2)
    sd = np.std(final_list2)
    final_list3=[x for x in final_list2 if (x >= mean - sd)]

    return np.e**(np.mean(final_list3))

def closestListInd(val, List):

    diffList=[l-val for l in List]
    absdifflist=[np.abs(l) for l in diffList]
    tmp=min(absdifflist)

    return absdifflist.index(tmp)

def getTprofile(yind,TTab,pars):

    Y=pars[0][0]
    """ input (horizontal index ) """
    if yind>=Y:
        print('error: yind must be strictly less than Y')

    return [float(el) for el in TTab[yind::Y]]
		
def getVprofile(yind,VTab,pars):

    Y=pars[0][0]
    """ input (horizontal index ) """
    if yind>=Y:
        print('error: yind must be strictly less than',Y)

    return [float(el) for el in VTab[yind::Y]]

def getInterpVal(val,X,Y,baseval):

    Y.insert(0,baseval)
    X.insert(0,0)

    return np.interp(val,X,Y)
    #interpval=1
    #for i in range(len(X)-1):
    #    if Y[i]>Y[i-1] and Y[i]>val and Y[i-1]<val:
    #        interpval= X[i-1]+((X[i]-X[i-1])/(Y[i]-Y[i-1]))*(val-Y[i-1])
    #    elif Y[i]<Y[i-1] and Y[i]<val and Y[i-1]>val:
    #        interpval = X[i]-(((X[i]-X[i-1])/(Y[i-1]-Y[i]))*(val-Y[i]))
    #if interpval==1:
    #    interpval = baseval+((X[0]-baseval)/(Y[0]-0))*(val-0)
    #return interpval

def getRadSlice(zind,Tab,pars):

    Y,Z=pars[0]
    start=zind*Y
    end=(zind+1)*Y
    #print('zind:',zind,'Y:',Y,' Z:',Z,' start:',start,' end:',end,'tab[0]: ',Tab[0], 'len tab:',len(Tab))
    if zind>=Z:
        print('error: zind must be strictly less than', Z)

    #print(len([float(el) for el in Tab[start:end]]))
    return [float(el) for el in Tab[start:end]]
	
def getRadSlice_ITPL(ri,ri_is_list,Tab,pars):

    Y,Z=pars[0]

    if ri_is_list:
        i0 = [int(np.floor(rr)) for rr in ri]
        i1 = [int(np.ceil(rr)) for rr in ri]

        start0=[b*Y for b in i0]
        end0=[(b+1)*Y for b in i0]
        start1=[b*Y for b in i1]
        end1=[(b+1)*Y for b in i1]

        val = []

        for horz_idx in range(len(ri)):

            thisr = ri[horz_idx]
            s0 = start0[horz_idx]
            e0 = end0[horz_idx]
            s1 = start1[horz_idx]
            e1 = end1[horz_idx]

            if thisr>=Z:
                print('error: zind must be strictly less than', Z)
            elif i0[horz_idx]==i1[horz_idx]:
                val.append(float(Tab[s0:e0][horz_idx]))
            else:
                x0 = float(Tab[s0:e0][horz_idx])
                x1 = float(Tab[s1:e1][horz_idx])
                val.append(x0+(thisr-i0[horz_idx])*(x1-x0)/(i1[horz_idx]-i0[horz_idx]))
        
        return val

    else:

        i0 = int(np.floor(ri))
        i1 = int(np.ceil(ri))
        start0=i0*Y
        end0=(i0+1)*Y
        start1=i1*Y
        end1=(i1+1)*Y

        if ri>=Z:
            print('error: zind must be strictly less than', Z)
        if i0==i1:
            val = [float(el) for el in Tab[start0:end0]]
        else:
            x0 = [float(el) for el in Tab[start0:end0]]
            x1 = [float(el) for el in Tab[start1:end1]]
            val = [x0[i]+(ri-i0)*(x1[i]-x0[i])/(i1-i0) for i in range(Y)]

        return val
		
def getTimeARR(direct,spt0,spt1,step):

    TIME=[]

    with open(direct+'_rprof.dat') as f_rprof:
        PROFreader = csv.reader(f_rprof)

        for row in PROFreader:
            if row!=[]:
                if row[0][0]=="*":
                    TIME.append(float(row[0][-10:])/31556952/10**6)

    return TIME[spt0:spt1+1:step]
		
# GET Z-values of BL THICKNESSESS INTERPOLATIONS
def getBLptsITPL(prof,rad,t,tcmb):

    flaghigh=0
    yhigh=0
    RI=0
    i=0

    while flaghigh==0:

        if prof[i]<t and i==0:
            flaghigh=1
            [y1,y2]= [0,rad[0]]
            [r1,r2]= [0,1]
            [x1,x2]= [tcmb,prof[0]]
            yhigh = (y2-y1)/(x2-x1)*(t-x1)+y1
            RI = (r2-r1)/(x2-x1)*(t-x1)+r1
            break

        if prof[i]<t:
            flaghigh=1
            [y1,y2]= [rad[i-1],rad[i]]
            [r1,r2]= [i-1,i]
            [x1,x2]= [prof[i-1],prof[i]]
            yhigh = (y2-y1)/(x2-x1)*(t-x1)+y1
            RI = (r2-r1)/(x2-x1)*(t-x1)+r1

        i=i+1

    return yhigh,RI

def getBLptsITPL_cold(prof,rad,t,tstag,convregthk,tcmb):

    flaglow=0
    ylow=0
    RI=0
    i=1

    while flaglow==0:

        if prof[i]<t:

            flaglow=1
            [y1,y2]= [rad[i-1],rad[i]]
            [r1,r2]= [i-1,i]
            [x1,x2]= [prof[i-1],prof[i]]

            ylow = convregthk-((y2-y1)/(x2-x1)*(t-x1)+y1)
            botstag=getInterpVal(tstag,np.arange(len(prof)),prof,tcmb)
            RI = botstag-((r2-r1)/(x2-x1)*(t-x1)+r1)

        i=i+1

    return ylow,RI
	
def getRAD(direct,pars):

    Z=pars[0][1]
    RAD=[]

    with open(direct+'_rprof.dat') as f_rprof:
        PROFreader = csv.reader(f_rprof, delimiter=" ", skipinitialspace=True)
        count=0

        for row in PROFreader:
            if row!=[]:
                if row[0]!='' and  row[0][0]!="*" and row[0][0]!='r' and count<Z:
                    RAD.append(float(row[0])/10**3)
                    count=count+1
    
    return RAD

def getRAlslice(lp,lkm,dt,tabV,pars):

    Y,Z=pars[0]
    a_exp,g,Cp,kappa,rho = pars[1][:5]
    i0=int(np.floor(lp))
    i1=i0+1
    bl=lkm*1000
    start0=i0*Y
    end0=(i0+1)*Y
    start1=i1*Y
    end1=(i1+1)*Y

    v0 = [float(el) for el in tabV[start0:end0]]
    v1 = [float(el) for el in tabV[start1:end1]]
    v = [v0[i]+(lp-i0)*(v1[i]-v0[i])/(i1-i0) for i in range(Y)]

    v = np.array(v)
    RA= [(a_exp*rho*g*(bl)**3*dt)/(kappa*v[i]) for i in range(Y)]
    return RA

def getRAlslice_2(Lp,Lkm,dt,tabV,pars):
    Y,Z=pars[0]
    a_exp,g,Cp,kappa,rho = pars[1][:5]
    RA=[]
    V=[]

    for i in range(Y):

        i0=int(np.floor(Lp[i]))
        i1=i0+1
        bl=Lkm[i]*1000
        start0=i0*Y
        end0=(i0+1)*Y
        start1=i1*Y
        end1=(i1+1)*Y

        v0 = float(tabV[start0:end0][i])
        v1 = float(tabV[start1:end1][i])

        v = v0+(Lp[i]-i0)*(v1-v0)/(i1-i0)

        V.append(v)
        RA.append((a_exp*rho*g*(bl)**3*dt)/(kappa*v))
    return RA
	
def plumeNo_Th(thld,tbl_ht,STRIP,factor,pars):
    Y=pars[0][0]
    core_THICK=pars[2][2]
    lll=len(STRIP)

    counter=0
    valcounter=0
    ind=0
    thks=[]
    inplume=0
    prev_val=1
    stepdist=(core_THICK+tbl_ht*1000*factor)*(np.pi/Y)/1000

    for val in STRIP:
        m=(val-prev_val)/stepdist

        if val!=prev_val:
            x=(thld-prev_val)/m
        else:
            x=0

        if inplume==0 and val>=thld and valcounter==0:
            counter=counter+1
            inplume=stepdist
        elif inplume==0 and val>=thld and valcounter!=0:
            counter=counter+1
            inplume=stepdist-x

        if inplume!=0 and val<thld:
            inplume=inplume+x
            thks.append(inplume)
            inplume=0
        if inplume!=0 and val>=thld:
            inplume=inplume+stepdist
        if ind==lll-1 and inplume!=0:
            thks.append(inplume+stepdist)
            inplume=0

        ind=ind+1
        valcounter=valcounter+1
        prev_val=val

    return [counter, thks]
	
def getSlabThick(centerTPROF,TGRAD,RAD,TLOW):
    inslab=0
    slabthick=0
    RAD=RAD[:200]

    for r in range(len(RAD)):
        if r==0 and centerTPROF[r]<TLOW:
            R0=0
            inslab=1

        elif centerTPROF[r]<TLOW and inslab==0:
            r0=r
            r0_prev=r-1
            Tprev=centerTPROF[r0_prev]
            T=centerTPROF[r0]
            m=(T-Tprev)/(RAD[r0]-RAD[r0_prev])
            R0 = (1/m)*(TLOW-T)+RAD[r0]
            inslab=1

        elif centerTPROF[r]>=TLOW and inslab==1:
            r1=r
            r1_prev=r-1
            Tprev=centerTPROF[r1_prev]
            T=centerTPROF[r1]
            m=(T-Tprev)/(RAD[r1]-RAD[r1_prev])
            R1 = (1/m)*(TLOW-T)+RAD[r1]
            slabthick=R1-R0
            break
    return slabthick


def clumpedSlice(binSLC,pl_num,p_wids,c_num,c_lst):
    p_labs=[]
    #print("K: ",c_num)

    for i in range(pl_num):
        #print('i: ',i,'wids: ',p_wids,'num:  ',pl_num,'c_num: ',c_num,'c_lst: ',c_lst)
        w = p_wids[i]
        l = 1000
        mindel=10000

        for j in range(c_num):

            if np.abs(c_lst[j]-w)<mindel:
                mindel=np.abs(c_lst[j]-w)
                l=j+1
        p_labs.append(l)

    #print("PLABELS: ",p_labs)
    #print("PWIDS: ",p_wids)
    outSLC=[]
    Y_HEAD_COORDS=[]
    ct_pls=0
    inpl=0
    ct=0

    for b in binSLC:

        if b==False and inpl==0:
            outSLC.append(0)
        elif b==True and inpl==0:
            #print(len(p_labs),ct_pls)
            outSLC.append(p_labs[ct_pls])
            if p_labs[ct_pls]==2:
                Y_HEAD_COORDS.append(ct)
            inpl=1
        elif b==True and inpl==1:
            outSLC.append(p_labs[ct_pls])
            if p_labs[ct_pls]==2:
                Y_HEAD_COORDS.append(ct)

        if b==False and inpl==1:
            outSLC.append(0)
            ct_pls=ct_pls+1
            inpl=0
        ct=ct+1

    return outSLC, Y_HEAD_COORDS

def getVrms(Tablevy,Tablevz,pars):
    Y,Z=pars[0]
    kappa=pars[1][3]

    VY_SQURD=np.array([(1/Y)*np.sum( np.array(getRadSlice(i,Tablevy))**2 ) for i in range(Z)])
    VZ_SQURD=np.array([(1/Y)*np.sum( np.array(getRadSlice(i,Tablevz))**2 ) for i in range(Z)])

    return np.sqrt(VY_SQURD+VZ_SQURD)

def getSTAGTHICK(tgrad,vgrad,rad,time,timen,modname,pars,tcmb,tt_str,showplotBool,showplotFreq):
    [Y,Z]=pars[0]
    [domain_THICK,air_THICK,core_THICK,D_MantlenStag,mantle_THICK]=pars[2]
    [a_exp,g,Cp,kappa,rho,therm_conductivity,H_flux,H_mass,ra_num_const]=pars[1]
    [mm0,mm1]=pars[3]

    zipT=list(zip(rad,tgrad))
    zipV=list(zip(rad,vgrad))
    zipT.reverse()
    zipV.reverse()

    stagthick=0

    t_mean = trimmedMean(tgrad[mm0:mm1])
    v_mean = trimmedMeanV(vgrad[mm0:mm1])
    t_stag = np.mean(tgrad[mm0:mm1])

    t_tbl=0

    v_limit=10*v_mean #np.min(vgrad[mm0:mm1])

    for d in range(0,Z,1):
        d1=d+1

        if v_limit>1e25:
            if zipV[d][1]>9.9e24 and zipV[d1][1]<9.9e24:
                r_stag=zipV[d][0]
                stagthick=(domain_THICK-air_THICK)/10**3-r_stag
                t_stag = (zipT[d][1] + tcmb)/2
                break

        if zipV[d][1]>v_limit and zipV[d1][1]<v_limit: #10*np.min(vgrad[mm0:mm1]):
            m = (v_limit-zipV[d][1])/(zipV[d1][1]-zipV[d][1])
            r_stag=zipV[d][0] - m*(zipV[d][0]-zipV[d1][0])
            stagthick=(domain_THICK-air_THICK)/10**3-r_stag
            t_stag = zipT[d][1] + m*(zipT[d1][1]-zipT[d][1])

            break
    
    Afont = {'fontname':'Arial','size':15}
    if showplotBool and timen%showplotFreq==0:

        fig, ax = plt.subplots(1)
        fig.set_size_inches(6,5)
        ax.set_title(str(time)[:5]+' Myrs. mantle temp:'+str(t_mean)[:5]+' \n thickness: '+ str(stagthick)[:5]+' temp @ bottom:'+str(t_stag)[:6])
        ax.plot(v_mean*np.ones([Z]),rad,c='teal',linewidth=3)
        ax.plot(vgrad,rad,c='k',linewidth=5)
        ax.plot(v_limit*np.ones([Z]),rad,c='goldenrod',linewidth=3,linestyle='-')
        ax.set_xlim([1e18,1e26])
        ax.set_xscale('log')
        ax.set_yticks(ticks=[0,500,1000,1500,2000,2500,3000],labels=[0,500,1000,1500,2000,2500,3000],fontsize=15)
        ax.set_xticks(ticks=[1e18,1e20,1e22,1e24,1e26],labels=[1e18,1e20,1e22,1e24,1e26],fontsize=15)
        ax.set_ylim([0,mantle_THICK])
        ax.grid(axis='both', which='major', color='0.25', linestyle='-')
        ax.grid(axis='x', which='minor', color='0.8', linestyle='--')
        plt.subplots_adjust(top=0.85)
        #plt.savefig('/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/PUB_FIG/v_prof/'+modname+'/'+str(timen)+'.png')
        plt.show()
        plt.clf()
        plt.close()

    
    convectingregionthickness = domain_THICK/10**3 - air_THICK/10**3 - stagthick  
    return convectingregionthickness, t_stag


def makePlumeScatterPlot(T,plumewids,plumeno):

    if len(T)!=len(plumeno) or len(T)!=len(plumewids):

        if len(T)-1==len(plumeno) and len(T)-1==len(plumewids):
            T=T[:-1]

        else:
            print('T:',len(T), 'no:', len(plumeno), 'thk', len(plumewids))
            raise NameError
    
    scatterpoints=[]
    ct=0

    for t in T:
        for n in range(plumeno[ct]):
            scatterpoints.append([t,plumewids[ct][n]])
        ct=ct+1
    
    t_plot=[pt[0] for pt in scatterpoints]
    w_plot=[pt[1] for pt in scatterpoints]

    fig, ax = plt.subplots(1)
    fig.set_size_inches(10,6)
    ax.scatter(t_plot,w_plot)
    ax.set_yscale('log')
    plt.show()
    plt.clf()
    plt.close()