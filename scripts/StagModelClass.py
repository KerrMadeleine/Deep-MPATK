# load necessary packages

from matplotlib.cbook import flatten
from sklearn.cluster import KMeans
from sklearn import metrics as skmet;
import numpy as np
import pandas as pd
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import csv
import os
import sys
import stagyypythonmodule as spm
from cmcrameri import cm


# Class structure to handle 2D StagYY models in spherical annulus

class StagYYModel:
    """
    highres_path
    longtime_path
    nmax
    nmin
    mesh_path
    filename
    modeldim_str
    modeldim_z
    modeldim_y
    Tcmb
    Tsurf
    TM_0
    output_path
    Rp
    Rc
    airthick
    alpha
    g
    Cp
    thermal_diffusivity
    ref_density
    thermal_conductivity
    cores_nz_0
    cores_ny_0
    cores_nz_f
    cores_ny_f
    cores_nzny_swap
    """

    def __init__(self, path, mesh_path, dim_y, dim_z, nmax, tcmb, Tm_0, tsurf, Rc, D):
        """
        initialize a StagYY model instance with a model path, mesh path, y and z dimensions
        maximum number of timesteps, CMB temperature, mantle temperature, surface temperature
        """
        self.path=path
        self.dim_y = dim_y
        self.dim_z = dim_z
        self.nmin=0
        self.nmax=nmax
        self.Tcmb=tcmb
        self.Tsurf=tsurf
        self.mesh_path= mesh_path + 'STAGmesh'+str(dim_y)+'x'+str(dim_z)+'.csv'
        self.Tm_0 = Tm_0
        self.Rc = Rc
        self.D = D 
        self.Rp = Rc + D

        self.r_prof_path = path+'_rprof.dat'
        self.TimeEvoDFs_dict = {}

        Tablemesh=[]
        with open(self.mesh_path) as f:
            reader = csv.reader(f)
            for col in reader:
                Tablemesh.append(col)
        self.y=[float(i) for i in Tablemesh[0][0:self.dim_y]] # get the radii coords
        self.z=[float(i) for i in Tablemesh[0][self.dim_y::2*self.dim_y]] # get the z coords
        self.Y=len(self.y)
        self.Z=len(self.z)
        self.arrY=[i for i in range(1,self.Y+1)]
        self.arrZ=[i for i in range(1,self.Z+1)]

        print('initializing ...')
        print('path:', self.path)
        print('dims:', self.dim_y, ' x ', self.dim_z)
        print('n = ', self.nmin, '- > ', self.nmax)
        print('T :', self.Tsurf, '-->', self.Tm_0, '-->', self.Tcmb)

        self.SEC_IN_YR = 31556926


    # ------------ StagYY Binary-------------------------------------


    # ---------------------------------------------------------------

    def setPhysPars(self, THERM_EXPANS, GRAV_ACCEL,SPECF_HEAT_CONST_PRESS,THERM_DIFFUS,REF_DENSITY,THERM_COND):
        self.alpha = THERM_EXPANS
        self.g = GRAV_ACCEL
        self.Cp = SPECF_HEAT_CONST_PRESS
        self.kappa = THERM_DIFFUS
        self.rho_0 = REF_DENSITY
        self.k_cond = THERM_DIFFUS*REF_DENSITY*SPECF_HEAT_CONST_PRESS

    def setOnsetPars(self,str_label):
        this_DF = self.TimeEvoDFs_dict[str_label]
        self.Ra_crit_hot = this_DF['Ra_crit_hot'][0]
        self.Ra_crit_cold = this_DF['Ra_crit_cold'][0]
        self.t_onset_hot = this_DF['time_onset_hot (Myrs)'][0]
        self.t_onset_cold = this_DF['time_onset_cold (Myrs)'][0]
        self.z_onset_hot = this_DF['hot_TBL_onset (m)'][0]
        self.z_onset_cold = this_DF['cold_TBL_onset (m)'][0]


        
    # ------------ StagYY tableT tableV------------------------------
    def getTableT(self,ts):
        """
        for a timestep ts, return an array nr x nphi of temperature values
        """
        ts = str(ts).zfill(5)
        openFileT= self.path+'t_'+str(ts)+'.csv' #temp file to open
        TableT=np.loadtxt(openFileT)
        return TableT
    
    def getTableV(self,ts):
        """
        for a timestep ts, return an array nr x nphi of temperature values
        """
        ts = str(ts).zfill(5)
        openFileV= self.path+'eta_'+str(ts)+'.csv' #visc file to open
        TableV=np.loadtxt(openFileV)
        return TableV

    


    # ---------------------------------------------------------------


    # ------------ Plotting Fields --------------------------------

    def plotT_2Dcart(self,T_table):
        """
        plot temperature in 2d cartesian
        """
        fig, ax = plt.subplots(1)
        fig.set_size_inches(9,3)
        c=ax.pcolor(self.arrY,self.arrZ,T_table-self.Tm_0,cmap=cm.vik, vmin=-200,vmax=200) #define color bar bounds with vmin and vmax
        cax = plt.axes([0.92, 0.25, 0.02, 0.5])
        fig.colorbar(c,cax=cax)
        plt.show()
        plt.clf()
        plt.close()
    
    def plotlogV_2Dcart(self,V_table):
        """
        plot log viscosity in 2d cartesian
        """
        fig, ax = plt.subplots(1)
        fig.set_size_inches(9,3)
        c=ax.pcolor(self.arrY,self.arrZ,np.log10(V_table),cmap=cm.oslo_r, vmin=18,vmax=25) #define color bar bounds with vmin and vmax
        cax = plt.axes([0.92, 0.25, 0.02, 0.5])
        fig.colorbar(c,cax=cax)
        plt.show()
        plt.clf()
        plt.close()

    # ---------------------------------------------------------------



    # -----------------------Radial profiles------------------------
    def getRadialprofile(self, table, y_idx):
        return table[:,y_idx]

    def getHorizontalprofile(self, table, z_idx):
        return table[z_idx]

    def getHorizontalprofile_z(self, table_str, ts, z_slice):
        """
                 _ _ _ _ _ _
                |_|_|_|_|_|_|
           > -- |_|_|_|_|_|_| -- <-
                |_|_|_|_|_|_|

        """

        if table_str == 'T':
            table = self.getTableT(ts)
        elif table_str == 'V':
            table = np.log10(self.getTableV(ts))
        else:
            raise ValueError('table_str not recognized')

        # print('SIZE:',self.Y, np.shape(table),self.path)

        z_onset_index_float = np.interp(z_slice, self.z, self.arrZ)
        z_onset_floor = int(np.floor(z_onset_index_float))
        z_onset_ceil = int(np.ceil(z_onset_index_float))

        slice_floor = table[z_onset_floor]
        slice_ceil = table[z_onset_ceil]

        # get the weighted mean between z floor and z ceil to get value at z_nset
        percent_floor_2_ceil = (z_onset_index_float - z_onset_floor)/(z_onset_ceil-z_onset_floor)
        slice_onset = slice_floor + percent_floor_2_ceil*(slice_ceil - slice_floor)

        # fig, ax = plt.subplots(1)
        # fig.set_size_inches(9.5,7)
        # plt.title(str(percent_floor_2_ceil))
        # plt.plot(self.arrY,slice_floor,color='r')
        # plt.plot(self.arrY,slice_ceil,color='b')
        # plt.plot(self.arrY,slice_onset,color='g')
        # plt.show()
        # plt.clf()
        # plt.close()

        return slice_onset

    def get_Ra_loc_hot_slice(self, ts, z_onset_hot):

        Tmean_evo = self.getmeanTevolution(2, 5)
        T_z_onset_hot = Tmean_evo[ts]#self.getHorizontalprofile_z('T', ts, z_onset_hot)
        V_z_onset_hot = 10**(self.getHorizontalprofile_z('V', ts, z_onset_hot))

        Ra_loc_star_slice = (self.rho_0 * self.g * self.alpha * (self.Tcmb - T_z_onset_hot) * z_onset_hot**3)/(self.kappa * V_z_onset_hot)
        
        # fig, ax = plt.subplots(1)
        # fig.set_size_inches(9.5,7)
        # plt.title(str(self.Ra_crit_hot))
        # plt.plot(self.arrY,Ra_loc_star_slice,color='k')
        # plt.plot(self.arrY,self.Ra_crit_hot*np.ones(self.Y),color='b')
        # plt.yscale('log')
        # plt.show()
        # plt.clf()
        # plt.close()

        return Ra_loc_star_slice


    def getMeanRadialprofile(self, table):
        for z_i in range(self.Z):
            print(self.getHorizontalprofile(table, z_i))
        profile = [np.mean(self.getHorizontalprofile(table, z_i)) for z_i in range(self.Z)]
        return profile

    def getMeanHorizontalprofile(self, table):
        profile = [np.mean(self.getRadialprofile(table, y_i)) for y_i in range(self.Y)]
        return profile

    def trimmedMeanSig(self,a_array):
        mu = np.mean(a_array)
        sigma = np.std(a_array)
        dev = a_array - mu
        trim = [a_array[i] for i in range(len(a_array)) if np.abs(dev[i])<sigma]
        if trim == []:
            return np.mean(a_array)
        else:
            return np.mean(trim)

        
    # def trimmedMean(self,a_array):
    #     mu = np.mean(a_array)
    #     sigma = 0.9 * (np.max(np.abs(mu-a_array)))
    #     dev = a_array - mu
    #     trim = [a_array[i] for i in range(len(a_array)) if np.abs(dev[i])<sigma]
    #     if trim == []:
    #         return np.mean(a_array)
    #     else:
    #         return np.mean(trim)

    def getTrimmedMeanHorizontalprofile(self, table):
        profile = [self.trimmedMean(self.getRadialprofile(table, y_i)) for y_i in range(self.Y)]
        return profile


    def getTrimmedMeanRadialprofile(self,table):
        profile = [self.trimmedMean(self.getHorizontalprofile(table, z_i)) for z_i in range(self.Z)]
        return profile

    #-----------------------------------------------------------------


    # ------------ Reading from rprof file. Profile plotting  ------

    def get_rprof_dictionary(self, par_n0, par_n1):
        """
        reads the _rprof.dat file
        outputs a dictionary with headers as keys and (time x depth) arrays as values
        """
        print('making dictionary from:', self.r_prof_path)

        TableData = list()
        TimeSepData = list()

        # read file into a list of lists
        with open(self.r_prof_path) as f:
            reader = csv.reader(f, delimiter=" ", skipinitialspace=True)
            headers = next(reader,None)
            #print('len headers:', len(headers))
            for row in reader:
                TableData.append(row)
        #print(headers[par_n0:par_n1])

        time = list()
        data_arr = list()
        this_time_data = list()

        #extract the time data and the measured data
        for row in TableData:
            if row[0][0]=='*': # if the row gives us a time   
                data_arr.append(this_time_data)             
                time.append(float(row[5])/self.SEC_IN_YR) #yrs
                this_time_data = list()
            else: # the row is not a time
                this_time_data.append([float(r) for r in row[par_n0:par_n1]])
        data_arr.append(this_time_data) #add the last set
    
        data_arr = data_arr[1:] # remove the empty list we appended in the first loop iteration
        data_arr = np.array(data_arr)
        time_arr = np.array([time] * len(data_arr[0]))
        time_arr = np.transpose(time_arr, (1,0))
        # transpose data to be like npar x ntime x ndepth
        data_arr = np.transpose(data_arr, (2,0,1))
        print(np.shape(data_arr), np.shape(time_arr))

        # make a dictionary based on the parameter names
        this_data_dict = {}
        this_data_dict['time'] = time_arr[:self.nmax+1] 
        for i_head in range(par_n0,par_n1):
            this_data_dict[headers[i_head]] = data_arr[i_head][:self.nmax+1] 

        return this_data_dict

    def set_rprof_dictionary(self, par_n0, par_n1): 
        """
        use get_rprof_dictionary to set the dictionary as a class attribute
        """
        self.rprof_data_dict = self.get_rprof_dictionary(par_n0,par_n1)

    def getTimeMyrs(self):
        """
        get a time evolution in Myrs
        """

        data_dict = self.rprof_data_dict

        time = data_dict['time'][:,0] * 1e-6 # yrs to myrs

        return time

    def getZkm(self):
        """
        get depths 1D array in km
        """

        data_dict = self.rprof_data_dict

        z = data_dict['r'][0] * 1e3 # m to km

        return z


    def getT_zon_ton(self, z_on,t_on_idx):
        """
        at the time index t_on_idx, what is the mean value of Temperature at the point z_on above the CMB
        """
        data_dict = self.rprof_data_dict
        T_profile = data_dict['Tmean'][t_on_idx]
        z_prof = data_dict['r'][t_on_idx]
        
        T_zon_ton = np.interp(z_on, z_prof,T_profile)

        # fig, ax = plt.subplots(1)
        # fig.set_size_inches(9.5,7)
        # plt.plot(T_profile, data_dict['r'][t_on_idx])
        # plt.plot(np.linspace(min(T_profile),max(T_profile), 10), np.ones(10) * z_on)
        # plt.show()
        # plt.clf()
        # plt.close()

        return T_zon_ton

    def getV_zon_ton(self, z_on,t_on_idx):
        """
        at the time index t_on_idx, what is the mean value of Viscosity at the point z_on above the CMB
        """
        data_dict = self.rprof_data_dict
        V_profile = data_dict['etalog'][t_on_idx]
        z_prof = data_dict['r'][t_on_idx]

        V_zon_ton = np.interp(z_on, z_prof,V_profile)

        print(self.path, self.Y, self.dim_y)
        V_slice = 10**(self.getHorizontalprofile_z('V', t_on_idx, z_on))

        # fig, ax = plt.subplots(1)
        # fig.set_size_inches(9.5,7)
        # plt.plot(V_profile, data_dict['r'][t_on_idx])
        # plt.plot(np.linspace(min(V_profile),max(V_profile), 10), np.ones(10) * z_on)
        # plt.xscale('log')
        # plt.show()
        # plt.clf()
        # plt.close()

        # fig, ax = plt.subplots(1)
        # fig.set_size_inches(9.5,7)
        # plt.plot(self.arrY, V_slice)
        # plt.plot(self.arrY, np.ones(self.Y) * V_zon_ton)
        # plt.yscale('log')
        # plt.show()
        # plt.clf()
        # plt.close()

        return V_zon_ton


    def make_T_profile_plot(self,ts,perc_dev,str_label):
        """
        makes a figure of the temperature profile at ts showing the cold and hot TBL thicknesses
        """
        data_dict = self.rprof_data_dict
        z = data_dict['r'][0]/1e3
        dep = self.D - z
        l = len(z)
        this_DF = self.TimeEvoDFs_dict[str_label]
        T_core = self.Tcmb
        T_SL = this_DF['T_SL_base (K)'][ts]
        T_mean = this_DF['Tmean (K)'][ts]

        T_cold_TBL = this_DF['cold_TBL_T (m)'][ts]
        T_hot_TBL = this_DF['hot_TBL_T (m)'][ts]

        z_cold_TBL = this_DF['cold_TBL_z (m)'][ts]/1e3
        z_hot_TBL = this_DF['hot_TBL_z (m)'][ts]/1e3

        SL_thick = this_DF['SL_thick (m)'][ts]
        SL_z = (self.D - SL_thick)/1e3


        fig, ax = plt.subplots(1)
        fig.set_size_inches(9.5,7)
        ax.set_title(str(this_DF['time (Myrs)'][ts]) + '  Myrs \n' + 'SL_thick = ' + str(SL_thick/1e3) + '  km')
        ax.plot([T_core] + list(data_dict['Tmean'][ts]), [0] + list(z), color='k',linewidth = 3)
        ax.plot(np.ones(10) * T_SL, np.linspace(z[0],z[-1],10), color='b',linewidth = 3)
        ax.plot(np.linspace(700,2000,10), np.ones(10) * SL_z, color ='b',linewidth = 3)
        ax.plot(np.linspace(700,2000,10), np.ones(10) * z_cold_TBL, color ='k',linewidth = 2,linestyle='--')
        ax.plot(np.linspace(700,2000,10), np.ones(10) * z_hot_TBL, color ='k',linewidth = 2,linestyle='--')
        ax.plot(np.ones(10) * T_core, np.linspace(z[0],z[-1],10), color='r',linewidth = 3)
        ax.plot(np.ones(10) * T_cold_TBL, np.linspace(z[0],z[-1],10), color='b',alpha=0.5,linewidth = 3)
        ax.plot(np.ones(10) * T_hot_TBL, np.linspace(z[0],z[-1],10), color='r',alpha=0.5,linewidth = 3)
        ax.set_xlim([690,2010])
        ax.set_ylim([0,2942])
        ax.set_ylabel('z above CMB (km)',fontsize=18)
        ax.set_xlabel('T',fontsize=18)
        ax.grid()
        plt.show()
        plt.clf()
        plt.close()


    def make_V_profile_plot(self,ts,eta_max,ord_of_mag_diff,str_label):

        """
        make a figure of the viscosity profile at ts showing the SL thickness
        """
        data_dict = self.rprof_data_dict
        z=data_dict['r'][0]/1e3 #km
        l = len(data_dict['r'][0])
        UM_min_eta = np.min(data_dict['etalog'][ts][l//2:])
        SL_eta = 10**ord_of_mag_diff * UM_min_eta
        z_eta_max = np.interp(.99*eta_max, data_dict['etalog'][ts][l//2:], z[l//2:])

        this_DF = self.TimeEvoDFs_dict[str_label]
        SL_thick = this_DF['SL_thick (m)'][ts]
        SL_z = (self.D - SL_thick)/1e3

        fig, ax = plt.subplots(1)
        fig.set_size_inches(5,7)
        ax.set_title(str(this_DF['time (Myrs)'][ts]) + '  Myrs \n' + 'SL_thick = ' + str(SL_thick/1e3) + '  km')
        ax.plot(data_dict['etalog'][ts], z, color = 'k',linewidth = 3)
        ax.plot(np.ones(10) * UM_min_eta, np.linspace(z[0],z[-1],10), color = 'gray',linewidth = 3)
        ax.plot(np.ones(10) * 10**ord_of_mag_diff  * UM_min_eta, np.linspace(z[0],z[-1],10), color='purple',linewidth = 3)
        ax.plot(np.linspace(1e18,1e26,10), np.ones(10) * SL_z,color='purple',linewidth = 3)
        ax.set_xscale('log')
        ax.set_xlim([9e17,1.2e25])
        ax.set_ylim([0,2942])
        ax.set_ylabel('z above CMB (km)')
        ax.set_xlabel('log[$\eta$ Pa•s] ')
        plt.grid()
        plt.show()
        plt.clf()
        plt.close()

    # -----------------------------------------------
    # -------- Evolution ----------------------------
    # -----------------------------------------------

    def getVrmsevolution(self, eights_0, eights_1):
        """
        get the rms velocity evolution in mm/yr within the range of eights_0*(L//8) - > eights_1*(L//8)
        where L is the number of points in the z direction (radial/depth).
        """
        # get the whole mantle VRMS
        data_dict = self.rprof_data_dict
        r = list(data_dict['r'][0]) + [self.D]
        areas = 0.5 * np.pi * (np.array(r) + self.Rc)** 2
        dA = np.diff(areas)

        l_8seg = len(r)//8

        r = data_dict['r'][0]
        r = r[l_8seg*eights_0:l_8seg*eights_1]
        dA = dA[l_8seg*eights_0:l_8seg*eights_1-1]
        domain_A =  0.5 * np.pi * ((max(r) + self.Rc)**2 - (min(r) + self.Rc)**2)
        VRMS = data_dict['vrms'] * self.SEC_IN_YR *1e3 # mm/yr
        wholedomain_vrms =   (1/domain_A) * np.array([sum(vrms_t[l_8seg*eights_0:l_8seg*eights_1-1]*dA) for vrms_t in VRMS])
     
        return wholedomain_vrms

    def getmeanTevolution(self, eights_0, eights_1):
        """
        get the mean temperature evolution in K within the range of eights_0*(L//8) - > eights_1*(L//8)
        where L is the number of points in the z direction (radial/depth).
        """
        data_dict = self.rprof_data_dict
        r = list(data_dict['r'][0]) + [self.D]
        areas = 0.5 * np.pi * (np.array(r) + self.Rc)** 2
        dA = np.diff(areas)
        l_8seg = len(r)//8
        r = data_dict['r'][0]
        r = r[l_8seg*eights_0:l_8seg*eights_1]
        dA = dA[l_8seg*eights_0:l_8seg*eights_1-1]
        domain_A =  0.5 * np.pi * ((max(r) + self.Rc)**2 - (min(r) + self.Rc)**2)

        Tmean = data_dict['Tmean'] 
        wholedomain_Tmean =   (1/domain_A) * np.array([sum(tmean_t[l_8seg*eights_0:l_8seg*eights_1-1]*dA) for tmean_t in Tmean])
     
        return wholedomain_Tmean

    def getmeanVevolution(self, eights_0, eights_1):
        """
        get the mean temperature evolution in K within the range of eights_0*(L//8) - > eights_1*(L//8)
        where L is the number of points in the z direction (radial/depth).
        """
        data_dict = self.rprof_data_dict
        r = list(data_dict['r'][0]) + [self.D]
        areas = 0.5 * np.pi * (np.array(r) + self.Rc)** 2
        dA = np.diff(areas)
        l_8seg = len(r)//8
        r = data_dict['r'][0]
        r = r[l_8seg*eights_0:l_8seg*eights_1]
        dA = dA[l_8seg*eights_0:l_8seg*eights_1-1]
        domain_A =  0.5 * np.pi * ((max(r) + self.Rc)**2 - (min(r) + self.Rc)**2)

        Vmean = np.log10(data_dict['etalog'])
        wholedomain_Vmean =  10 ** ((1/domain_A) * np.array([sum(vmean_t[l_8seg*eights_0:l_8seg*eights_1-1]*dA) for vmean_t in Vmean]))
     
        return wholedomain_Vmean 


    def getSLthickevolution(self, eta_max, ord_of_mag_diff):
        """
        return the stagnant lid thickness and the convecting region domain (D- SL_thick)
        evoltion over time

        the min viscosity of the SL is ord_of_mag_diff greater than the minimum radial average viscosity
        in the upper half of the domain.
        """
        data_dict = self.rprof_data_dict

        z = data_dict['r'][0]
        dep = self.D - z
        l = len(z)

        SL_thick_evo = list()
        Conv_D_evo = list()
        T_SL_base_evo = list()

        for t_idx in range(0,self.nmax+1):

            UM_min_eta = np.min(data_dict['etalog'][t_idx][l//2:])
            SL_eta = 10**ord_of_mag_diff * UM_min_eta
            z_eta_max = np.interp(.99*eta_max, data_dict['etalog'][t_idx][l//2:], data_dict['r'][t_idx][l//2:])

            if SL_eta > eta_max:
                SL_z = z_eta_max
                SL_thick = self.D - z_eta_max
            else:
                SL_z = np.interp(SL_eta, data_dict['etalog'][t_idx][l//2:], data_dict['r'][t_idx][l//2:])
                SL_thick = self.D - SL_z

            T_SL_base = np.interp(SL_thick, dep[::-1], data_dict['Tmean'][t_idx][::-1])
            
            # --- uncomment to see viscosity profiles ----------
            # fig, ax = plt.subplots(1)
            # fig.set_size_inches(9.5,7)
            # plt.plot(data_dict['etalog'][t_idx], data_dict['r'][t_idx])
            # plt.plot(np.ones(10) * np.min(data_dict['etalog'][t_idx][l//2:]), np.linspace(z[0],z[-1],10))
            # plt.plot(np.ones(10) * 10**ord_of_mag_diff  * np.min(data_dict['etalog'][t_idx][l//2:]), np.linspace(z[0],z[-1],10))
            # plt.plot(np.linspace(1e18,1e26,10), np.ones(10) * np.interp(.99*eta_max, data_dict['etalog'][t_idx][l//2:], data_dict['r'][t_idx][l//2:]))
            # plt.plot(np.linspace(1e18,1e26,10), np.ones(10) * SL_z, color ='k')
            # plt.xscale('log')
            # plt.show()
            # plt.clf()
            # plt.close()

            # --- uncomment to see temperature profiles ----------
            # fig, ax = plt.subplots(1)
            # fig.set_size_inches(9.5,7)
            # plt.plot(data_dict['Tmean'][t_idx], data_dict['r'][t_idx])
            # plt.plot(np.ones(10) * T_SL_base, np.linspace(z[0],z[-1],10))
            # plt.plot(np.linspace(700,2000,10), np.ones(10) * SL_z, color ='k')
            # plt.show()
            # plt.clf()
            # plt.close()

            SL_thick_evo.append(SL_thick)
            Conv_D_evo.append(SL_z)
            T_SL_base_evo.append(T_SL_base)

        return np.array(SL_thick_evo), np.array(Conv_D_evo), np.array(T_SL_base_evo)

    def getTimeOverturns(self,str_label):

        this_DF = self.TimeEvoDFs_dict[str_label]

        this_time = this_DF['time (Myrs)'] * 1e6 #yrs
        this_vrms = this_DF['vrms (mm/yr)'] * 1e-3 #m/yr
        this_D_conv = this_DF['D_conv (m)']     #m
        this_tau = this_D_conv/this_vrms

        tau = np.mean(this_D_conv[self.nmax//2:]/this_vrms[self.nmax//2:])

        return this_time/tau



    def getTBL_z_evolution(self, perc_dev, str_label):
        """
        takes in perc_dev to compute the percent deviation from mean T defines the TBL
        returns the distance from the core (z) from the hot TBL and cold TBL evolutions
        returns the temperature along the hot and cold TBLs
        returns the viscsoity along the hot and cold TBLs
        """

        data_dict = self.rprof_data_dict
        z = data_dict['r'][0]
        dep = self.D - z
        l = len(z)
        this_DF = self.TimeEvoDFs_dict[str_label]
        T_core = self.Tcmb
        T_SL = this_DF['T_SL_base (K)']
        T_mean = this_DF['Tmean (K)']
        V_mean = this_DF['Vmean (Pa s)']

        hot_TBL_z_evo = list()
        cold_TBL_z_evo = list()
        hot_TBL_T_evo = list()
        hot_TBL_V_evo = list()
        cold_TBL_T_evo = list()
        cold_TBL_V_evo = list()

        for t_idx in range(0,self.nmax+1):

            T_cold_TBL = T_SL[t_idx] + (1-perc_dev)*(T_mean[t_idx] - T_SL[t_idx])
            T_hot_TBL = T_mean[t_idx] + (perc_dev*(T_core - T_mean[t_idx]))

            z_hot_TBL = np.interp(-1*T_hot_TBL, [-1*T_core] + list(-1 * data_dict['Tmean'][t_idx][:l//2]), [0] + list(z[:l//2]))
            z_cold_TBL = np.interp(T_cold_TBL, data_dict['Tmean'][t_idx][::-1], z[::-1])

            V_cold_TBL = np.interp(z_cold_TBL, z[l//2:] , data_dict['etalog'][t_idx][l//2:])
            V_hot_TBL = np.interp(z_hot_TBL, z[:l//2] , data_dict['etalog'][t_idx][:l//2])

            hot_TBL_T_evo.append(T_hot_TBL)
            cold_TBL_T_evo.append(T_cold_TBL)
            hot_TBL_z_evo.append(z_hot_TBL)
            cold_TBL_z_evo.append(z_cold_TBL)
            hot_TBL_V_evo.append(V_hot_TBL)
            cold_TBL_V_evo.append(V_cold_TBL)

            # ----- plot temperature profiles ------------
            # fig, ax = plt.subplots(1)
            # fig.set_size_inches(9.5,7)
            # plt.plot(data_dict['Tmean'][t_idx] ,z)
            # plt.plot(T_SL[t_idx]*np.ones(10),np.linspace(z[0], z[-1], 10))
            # plt.plot(T_mean[t_idx]*np.ones(10),np.linspace(z[0], z[-1], 10))
            # plt.plot(T_core*np.ones(10),np.linspace(z[0], z[-1], 10))
            # plt.plot(T_cold_TBL*np.ones(10),np.linspace(z[0], z[-1], 10),color='k')
            # plt.plot(np.linspace(700,2000,10), np.ones(10) * z_cold_TBL, color='blue')
            # plt.plot(T_hot_TBL*np.ones(10),np.linspace(z[0], z[-1], 10), color='k')
            # plt.plot(np.linspace(700,2000,10), np.ones(10) * z_hot_TBL, color='red')
            # plt.show()
            # plt.clf()
            # plt.close()

            # ----- plot viscosity profiles ------------
            # fig, ax = plt.subplots(1)
            # fig.set_size_inches(9.5,7)
            # plt.plot(data_dict['etalog'][t_idx] ,z)
            # plt.plot(V_cold_TBL*np.ones(10), np.linspace(z[0],z[-1],10),color='blue')
            # plt.plot(np.linspace(1e18,1e25,10), np.ones(10) * z_cold_TBL, color='blue')
            # plt.plot(V_mean[t_idx]*np.ones(10),np.linspace(z[0], z[-1], 10),color='k')
            # plt.plot(V_hot_TBL*np.ones(10), np.linspace(z[0],z[-1],10),color='red')
            # plt.plot(np.linspace(1e18,1e25,10), np.ones(10) * z_hot_TBL, color='red')
            # plt.xscale('log')
            # plt.show()
            # plt.clf()
            # plt.close()
        
        return np.array(hot_TBL_z_evo), np.array(cold_TBL_z_evo), np.array(hot_TBL_T_evo), np.array(hot_TBL_V_evo), np.array(cold_TBL_T_evo), np.array(cold_TBL_V_evo)



    # -----------------------------------------------
    # -----------------------------------------------


    # ---------------------------------------------------------------------------

    def createNewEmptyTimeEvoDF(self,str_label):
        """
        create a new empty dataframe and associated csv file with the name '[str_label].csv'
        """

        if type(str_label) != str:
            raise ValueError('label needs to be a string')

        TimeEvoDF =  pd.DataFrame()

        filename = self.path + str_label + '.csv'

        print('check:', filename)

        if os.path.isfile(filename):
            raise ValueError(str_label+'.csv already exists')

        try:
            print(self.TimeEvoDFs_dict[str_label], ' already exists: str_label:', str_label)
        
        except:
            self.TimeEvoDFs_dict[str_label] = TimeEvoDF

        print('create:', filename)
        self.TimeEvoDFs_dict[str_label].to_csv(self.path+str_label+'.csv', index=False)



    def writeTimeEvoDF(self,str_label):
        """
        take a dataframe and write it to the csv file with the name '[str_label].csv'
        this is like "update csv with DataFrame changes"
        """

        filename = self.path + str_label + '.csv'
        print('write:', filename)
        try:
            self.TimeEvoDFs_dict[str_label].to_csv(self.path + str_label + '.csv', index=False)
        except:
            if not os.path.isfile(filename):
                self.createNewEmptyTimeEvoDF(str_label)
            else:
                self.TimeEvoDFs_dict[str_label] = pd.read_csv(self.path + str_label + '.csv')
                self.TimeEvoDFs_dict[str_label].to_csv(self.path + str_label + '.csv', index=False)

      
    def addSeriestoTimeEvoDF(self,str_label,series_label,series_data):
        """
        add a time series to a dataframe associated with the csv file '[str_label].csv'
        """
        filename = self.path + str_label + '.csv'
        if not os.path.isfile(filename):
            raise ValueError(str_label+'.csv does not exist')
        else:
            print(str_label)
            print(self.TimeEvoDFs_dict.keys())
            print(self.TimeEvoDFs_dict)
            this_time_evo_df = self.TimeEvoDFs_dict[str_label]
            prev_df_len = len(this_time_evo_df)
            new_df_len = len(series_data)

            if new_df_len>prev_df_len:
                this_time_evo_df = this_time_evo_df.reindex(range(new_df_len))
            elif prev_df_len > new_df_len:
                raise ValueError('Series too short!, prev:', prev_df_len, ' new:',new_df_len )

            this_time_evo_df[series_label]=series_data
            self.TimeEvoDFs_dict[str_label] = this_time_evo_df 
            print('now write')
            self.writeTimeEvoDF(str_label)
        

    def setTimeEvoDF(self, str_label):
        """
        initialize the dataframe object as a part of a class attribute (a dictionary of dataframes)
        """
        filename = self.path + str_label + '.csv'
        if not os.path.isfile(filename):
            raise ValueError(filename + 'does not exist')
        else:
            try:
                this_df = pd.read_csv(self.path + str_label + '.csv')
                print('Data loaded successfully')
            except:
                this_df = pd.DataFrame()
            self.TimeEvoDFs_dict[str_label] = this_df

    def readTimeEvoDF(self, str_label):
        """
        read in the csv file as a dataframe and returns this dataframe.
        """
        filename = self.path + str_label + '.csv'
        if not os.path.isfile(filename):
            raise ValueError(filename + 'does not exist')
        else:
            this_df = pd.read_csv(self.path + str_label + '.csv')
            return this_df
    # ---------------------------------------------------------------


    # ----------------- Extract from Data Frames ------------------

    def getRa_bw_t0_t1(self,str_label,ts_0,ts_1):
        
        this_DF = self.TimeEvoDFs_dict[str_label]

        Ra_WM_evo = this_DF['Ra_WM']

        return np.mean(Ra_WM_evo[ts_0:ts_1+1])

        

    # --------------------------------------------------------------



    # ------------ Detecting plume material ------------------------

    def plumeNo_Th(self,thld,tbl_ht,STRIP,factor):
        lll=len(STRIP)

        counter=0
        valcounter=0
        ind=0
        thks=[]
        inplume=0
        prev_val=1
        stepdist=(self.Rc+tbl_ht*factor)*(np.pi/self.Y)


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


        # # for debugging: 
        # fig, ax = plt.subplots(1)
        # fig.set_size_inches(9,4)
        # ax.set_title(str(counter)+'\n'+str(thks),fontsize=5)
        # ax.plot(self.y, STRIP)
        # ax.plot(self.y, thld * np.ones(self.Y))
        # plt.show()
        # plt.clf()
        # plt.close()

        return [counter, thks]

    # ---------------------------------------------------------------





    # ------------ Kellogg plots -----------------------------------

    def getColorMaps(self):
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

        return cmap_rb , cmap_bin


    def getPowerLawCoeffs(self, x_arr,y_arr):
        #if x_arr ==[] or y_arr==[]:
        #    return 0,0,[0]
        log_x = np.log(x_arr)
        log_y = np.log(y_arr)

        coeff, residuals, rank, singular_values, rcond = np.polyfit(x=log_x,y=log_y,deg=1,full=True)
        coeff, V = np.polyfit(x=log_x,y=log_y,deg=1,full=False,cov = True)
        print(np.polyfit(x=log_x,y=log_y,deg=1,full=True))
        print(np.polyfit(x=log_x,y=log_y,deg=1,full=False,cov = True))
        SS_res=residuals
        SS_tot = np.sum((log_y-np.mean(log_y))**2)
        Rsquared = 1 - (SS_res/SS_tot)

        return coeff[0],coeff[1],Rsquared, V

    def makeRaloc_star_Kellog(self):
        data_dict = self.rprof_data_dict
        time = data_dict['time'][:,0]
        Kellog_arr = list()
        for ts in range(self.nmax+1):
            print(ts,'/',self.nmax, '                    ', end='\r')
            Ra_slice = self.get_Ra_loc_hot_slice( ts, self.z_onset_hot)
            Kellog_arr.append(Ra_slice)
        Kellog_arr = np.array(Kellog_arr)

        #print(np.shape(self.arrY),np.shape(time), np.shape(Kellog_arr))

        fig, ax = plt.subplots(1)
        fig.set_size_inches(9,3)
        c=ax.pcolor(self.arrY,time,Kellog_arr,cmap=cm.devon_r, vmin=0,vmax=4000) #define color bar bounds with vmin and vmax
        cax = plt.axes([0.92, 0.25, 0.02, 0.5])
        fig.colorbar(c,cax=cax)
        plt.show()
        plt.clf()
        plt.close()

    def makeRaloc_star_Kellog_masked(self, Ra_crit):
        """
        Use a value of Ra_crit and the get_Ra_loc_hot_slice function
        to make a masked Montague-Kellog plot
        """
        data_dict = self.rprof_data_dict

        time = data_dict['time'][:,0]
        Kellog_arr = list()
        Mask_arr = list()
        for ts in range(self.nmax+1):
            print(ts,'/',self.nmax, '                    ', end='\r')
            Ra_slice = self.get_Ra_loc_hot_slice( ts, self.z_onset_hot)
            Kellog_arr.append(Ra_slice)
            Ra_bin_slice = np.array([x>=Ra_crit for x in Ra_slice])
            Mask_arr.append(Ra_bin_slice)
        Kellog_arr = np.array(Kellog_arr)
        Mask_arr = np.array(Mask_arr)

        #print(np.shape(self.arrY),np.shape(time), np.shape(Kellog_arr))

        cmap_rb , cmap_bin = self.getColorMaps()

        fig, ax = plt.subplots(1)
        fig.set_size_inches(9,3)
        c=ax.pcolor(self.arrY,time,Kellog_arr,cmap=cm.devon_r, vmin=0,vmax=4000) #define color bar bounds with vmin and vmax
        ax.pcolor(self.arrY,time,Mask_arr,cmap=cmap_bin, vmin=0,vmax=1)
        cax = plt.axes([0.92, 0.25, 0.02, 0.5])
        fig.colorbar(c,cax=cax)
        plt.show()
        plt.clf()
        plt.close()

    def saveRaloc_star_Kellog(self, Ra_crit, saveMaskbool, showplotbool):

        print('saving Ra_loc_* to text.....')
        data_dict = self.rprof_data_dict
        time = data_dict['time'][:,0]
        Kellog_arr = list()
        Mask_arr = list()
        for ts in range(self.nmax+1):
            print(ts,'/',self.nmax, '                    ', end='\r')
            Ra_slice = self.get_Ra_loc_hot_slice( ts, self.z_onset_hot)
            Kellog_arr.append(Ra_slice)
            Ra_bin_slice = np.array([x>=Ra_crit for x in Ra_slice])
            Mask_arr.append(Ra_bin_slice)
        Kellog_arr = np.array(Kellog_arr)
        Mask_arr = np.array(Mask_arr)

        np.savetxt(self.path+'Kellog_Ralocstar.txt',Kellog_arr)
        if saveMaskbool:
            np.savetxt(self.path+'Supercritical_Mask.txt',Mask_arr)
        
        if showplotbool:
            cmap_rb , cmap_bin = self.getColorMaps()
            fig, ax = plt.subplots(1)
            fig.set_size_inches(9,3)
            c=ax.pcolor(self.arrY,time,Kellog_arr,cmap=cm.devon_r, vmin=0,vmax=4000) #define color bar bounds with vmin and vmax
            if saveMaskbool:
                ax.pcolor(self.arrY,time,Mask_arr,cmap=cmap_bin, vmin=0,vmax=1)
            cax = plt.axes([0.92, 0.25, 0.02, 0.5])
            fig.colorbar(c,cax=cax)
            plt.show()
            plt.clf()
            plt.close()

    def readKellog_Ralocstar(self, plotBool):

        Ra_loc_star_Kellog = np.loadtxt(self.path+'Kellog_Ralocstar.txt')

        data_dict = self.rprof_data_dict
        time = data_dict['time'][:,0]

        if plotBool:
            fig, ax = plt.subplots(1)
            fig.set_size_inches(9,3)
            c=ax.pcolor(self.arrY,time,Ra_loc_star_Kellog,cmap=cm.devon_r, vmin=0,vmax=4000) #define color bar bounds with vmin and vmax
            #ax.pcolor(self.arrY,time,Mask_arr,cmap=cmap_bin, vmin=0,vmax=1)
            cax = plt.axes([0.92, 0.25, 0.02, 0.5])
            fig.colorbar(c,cax=cax)
            plt.show()
            plt.clf()
            plt.close()

        return Ra_loc_star_Kellog

    def readKellog_Mask(self,plotBool):

        Mask_Kellog = np.loadtxt(self.path+'Supercritical_Mask.txt')
        
        data_dict = self.rprof_data_dict
        time = data_dict['time'][:,0]
        
        if plotBool:
            cmap_rb , cmap_bin = self.getColorMaps()
            fig, ax = plt.subplots(1)
            fig.set_size_inches(9,3)
            #c=ax.pcolor(self.arrY,time,Ra_loc_star_Kellog,cmap=cm.devon_r, vmin=0,vmax=4000) #define color bar bounds with vmin and vmax
            c = ax.pcolor(self.arrY,time,Mask_arr,cmap=cmap_bin, vmin=0,vmax=1)
            cax = plt.axes([0.92, 0.25, 0.02, 0.5])
            fig.colorbar(c,cax=cax)
            plt.show()
            plt.clf()
            plt.close()

        return Mask_Kellog

    # ---------------------------------------------------------------







    def getPlumeCounts(self,ts_0,ts_1, Ra_loc_star_Kellog,showMiniplotsBool):
        all_plumes_snippet = list()
        number_of_all_plumes_snippet = 0

        each_plume_slice = list()
        number_of_each_plumes_slice = list()

        for ts in range(ts_0,ts_1):

            #print('time step = ', ts)

            plume_no_slice, plume_width_slice = self.plumeNo_Th(self.Ra_crit_hot ,self.z_onset_hot,Ra_loc_star_Kellog[ts],1)

            each_plume_slice.append(plume_width_slice)
            number_of_each_plumes_slice.append(plume_no_slice)

            all_plumes_snippet.extend(plume_width_slice)
            number_of_all_plumes_snippet = number_of_all_plumes_snippet + plume_no_slice 


        all_plumes_snippet=sorted(all_plumes_snippet)
        all_plumes_snippet=np.array(all_plumes_snippet)
        all_plumes_snippet.reshape(-1,1)

        return number_of_all_plumes_snippet, all_plumes_snippet, number_of_each_plumes_slice, each_plume_slice


    def Kmeans_all_plumes(self, kmax , allplumes):

        """
        takes in a 1D array of all plume widths (allplumes) and an ideal max number of clusters
        returns the centroids of each k and the cluster statistics
        """

        #  * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        # * * * * * * * *     KMEANS      * * * * * * * * * * * * 
        # * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

        kmax = np.min([kmax,len(allplumes)-1])


        CENTROIDS=[]
        CLUSTER_STDEVS=[]
        CLUSTER_NUMS=[]
        HD_OV_CON_SZ = list()
        HD_OV_CON_NO = list()
        AVG_PL_NO = list()

        data = allplumes.reshape(-1,1)
        stats_empty=1

        print('number of plumes:', len(allplumes), 'kmax = ', kmax)
        
        if len(allplumes)==0:
            return [[0, []]] , np.array([]),np.array([]),np.array([]),np.array([])
        else:
            SSdiff_CC=[]
            silo_1mod=[]
            DB_1mod=[]
            CH_1mod=[]
            GAP_1mod=[]

            for ii in range(1,kmax+1):
                print(ii,'/', kmax, end='                    \r')

                kmeans = KMeans(n_clusters=ii, random_state=0).fit(data)
                ssdiff=kmeans.inertia_
                labels=kmeans.fit_predict(data)

                std_per_feat=[]
                num_per_feat=[]

                for jj in range(ii):
                    std_per_feat.append(np.std([allplumes[x] for x in range(len(allplumes)) if labels[x]==jj]))
                    num_per_feat.append(len([allplumes[x] for x in range(len(allplumes)) if labels[x]==jj]))
                

                SSdiff_CC.append(ssdiff)

                if ii!=1:
                    # # DEBUG
                    # print('k = ',ii )
                    # print('data=',len(data), data)
                    # print('labels=',len(labels), labels)
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

                ## DEBUG
                #print(ii,":  ",means, std_per_feat, num_per_feat)

                clumpdata=list(zip(means,std_per_feat,num_per_feat))
                clumpdata=sorted(clumpdata)

                CENTROIDS.append([ii,[x[0] for x in clumpdata]])
                CLUSTER_STDEVS.append([ii,[x[1] for x in clumpdata]])
                CLUSTER_NUMS.append([ii,[x[2] for x in clumpdata]])

                #print('[',ii,',',[x[0] for x in clumpdata],',',[x[1] for x in clumpdata],',',[x[2] for x in clumpdata],']')

                # if ii>1:

                #     fig, ax = plt.subplots(1)
                #     fig.set_size_inches(5,3)
                #     fig.suptitle(str(int(self.Tm_0)))
                #     ax.scatter([x[0] for x in clumpdata],[x[2] for x in clumpdata])
                #     ax.errorbar([x[0] for x in clumpdata],[x[2] for x in clumpdata],xerr=[x[1] for x in clumpdata])
                #     plt.show()
                #     plt.clf()
                #     plt.close()

                if origDisp==0.0:
                    break

            # # DEBUG
            # print('Centroids:')
            # print(CENTROIDS)
            # print('statistics:')
            # print(SSdiff_CC, silo_1mod, DB_1mod, CH_1mod, GAP_1mod)

            return CENTROIDS, CLUSTER_STDEVS, CLUSTER_NUMS, np.array(SSdiff_CC), \
                np.array(silo_1mod) , np.array(DB_1mod), np.array(CH_1mod), \
                np.array(GAP_1mod)


    def clumpedSlice(self,binSLC,pl_num,p_wids,c_num,c_lst):
        p_labs=[]
        #print("total k: ",c_num)
        for i in range(pl_num):
            #print('i:',i,' wids:',p_wids,' num:',pl_num,' c_num:',c_num,'  c_lst:',c_lst)
            w = p_wids[i]
            l = 1
            mindel=10000000000

            #print('w = ', w)
            for j in range(c_num):
                #print('centroid number',j, ':',c_lst[j], ' mindel=', mindel, '  and l=',l)
                if np.abs(c_lst[j]-w)<mindel:
                    mindel=np.abs(c_lst[j]-w)
                    l=j+1
                    #print('diff = ', np.abs(c_lst[j]-w), 'and mindel=', mindel, '  and l=',l)
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

        # fig, ax = plt.subplots(1)
        # fig.set_size_inches(9,3)
        # ax.plot(self.arrY, outSLC)
        # plt.show()
        # plt.clf()
        # plt.close()

        return outSLC, Y_HEAD_COORDS


    def makeStatsPlots(self, kmeans_out):
        #print(kmeans_out)
        if kmeans_out != []:
            #print('kmeans_out:', kmeans_out[0], '\n   \n')
            (CENTROIDS_1mod,STDEV_1mod,NUM_1mod,SSdiff_CC, silo_1mod, DB_1mod, CH_1mod, GAP_1mod) = kmeans_out
            fig, ax = plt.subplots(6)
            fig.set_size_inches(5,7)
            fig.suptitle(str(int(self.Tm_0)))
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
            plt.show()
            plt.clf()
            plt.close()
        else:
            print('No plumes found at ts = ', ts, '\n ** \n')



    def plotGroupedKellog(self,ts_0,ts_1,k_choice,Ra_loc_star_Mask,number_of_each_plumes_slice,each_plume_slice,Centroids):
        
        k_cluster = len(Centroids)
        for cc in range(0,k_cluster,1):

            if k_choice>k_cluster:
                raise ValueError('bad cluster choice. k_choice must be < len(centroids)')

            if cc != k_choice -1:
                continue

            Ra_loc_star = list()
            Plume_groups = list()
            ts_list = list()
            for ts in range(ts_0,ts_1):
                labelled_slice, extrathing = self.clumpedSlice(Ra_loc_star_Mask[ts],number_of_each_plumes_slice[ts-ts_0],each_plume_slice[ts-ts_0],Centroids[cc][0],Centroids[cc][1])
                Ra_loc_star.append(Ra_loc_star_Mask[ts])
                Plume_groups.append(labelled_slice)
                ts_list.append(ts)


            cmap_rb , cmap_bin = self.getColorMaps()
            fig, ax = plt.subplots(1) # JAN4
            fig.set_size_inches(16,6)
            XCMB,TIMEevo = np.meshgrid(self.arrY,ts_list)
            ax.pcolor(XCMB,TIMEevo,Ra_loc_star,cmap='Blues')#,vmin=0,vmax=4000) #SPLIT
            c  = ax.pcolor(XCMB,TIMEevo,Plume_groups,cmap=cmap_rb)#,vmin=1,vmax=cc+1) # groups #SPLIT
            ax.set_ylim([ts_list[0]-1,ts_list[-1]+1])
            cax = plt.axes([0.93, 0.1, 0.02, 0.8])
            plt.show()
            plt.clf()
            plt.close()
            




    # ------------ Snippets     -----------------------------------


    def getHeadConduitStats(self,kmeans_out):
        Centroids = kmeans_out[0]
        Stdevs = kmeans_out[1]
        PlumesperCentroid = kmeans_out[2]

        k2_Cent = Centroids[1][1]
        k2_Stdevs = Stdevs[1][1]
        k2_ppC = PlumesperCentroid[1][1]

        print(k2_Cent, k2_Stdevs, k2_ppC)

        return (k2_Cent[0], k2_Stdevs[0], k2_ppC[0]), (k2_Cent[1], k2_Stdevs[1], k2_ppC[1])


    def getListOfAllHeads_n_Conds(self,allplumes):

        """
        takes in a 1D array of all plume widths (allplumes) and an ideal max number of clusters
        returns the centroids of each k and the cluster statistics
        """

        data = allplumes.reshape(-1,1)
        stats_empty=1
        
        if len(allplumes)==0:
            return [], []
        else:

            kmeans = KMeans(n_clusters=2, random_state=0).fit(data)
            ssdiff=kmeans.inertia_
            labels=kmeans.fit_predict(data)

            std_per_feat=[]
            num_per_feat=[]

            conduits_list = [allplumes[x] for x in range(len(allplumes)) if labels[x]==0]
            heads_list = [allplumes[x] for x in range(len(allplumes)) if labels[x]==1]

            return heads_list, conduits_list

    def create_snippet_bounds(self, ts_per_snippet, ts_start):

        endpoints = np.arange(ts_start,self.nmax, ts_per_snippet)

        bounds = [(endpoints[i],endpoints[i+1]) for i in range(len(endpoints)-1)]

        return bounds
    # ---------------------------------------------------------------