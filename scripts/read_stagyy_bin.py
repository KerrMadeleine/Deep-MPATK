import set_DeepMPATK_pars as par
import numpy as np
import matplotlib.pyplot as plt
import csv

stag_mesh_folder=par.STAGYY_OUTPUT_MESH_FOLDER
stag_out_folder=par.STAGYY_OUTPUT_FOLDER

for mod in par.STAGYY_OUTPUT_MODS[par.STAGYY_OUTPUT_MODS_0:par.STAGYY_OUTPUT_MODS_1]:
    print('converting model from bin -> csv for:', mod)
    ###########################################################################################
    # USER INPUTS:
    dir = stag_out_folder+mod  #ext_3500_1300_dp0_0/'
    fields = par.READSTAGBIN_fields # field name 't', 'eta', doesnt work with velocity/pressure yet.
    tstep0= par.READSTAGBIN_tstep0
    tstep1=par.READSTAGBIN_tstep1 # time-step
    tstep_step=par.READSTAGBIN_tstep_step
    y_dim = par.READSTAGBIN_y_dim # grid dim on horz axis
    z_dim = par.READSTAGBIN_z_dim # grid dim on vert axis
    nz_cores = par.READSTAGBIN_nz_cores #2 #number of core-blocks vertically
    ny_cores = par.READSTAGBIN_ny_cores #8 # number of core-blocks horizontally
    # path to csv created from staglab of the grid used in the run
    meshdir=stag_mesh_folder+'STAGmesh'+str(y_dim)+'x'+str(z_dim)+'.csv'
    showPlot_bool=par.READSTAGBIN_showPlot_bool
    writetocsv_bool=par.READSTAGBIN_writetocsv_bool
    ##########################################################################################

    nz_dim = z_dim//nz_cores #64
    ny_dim = y_dim//ny_cores #32

    for field in fields: # for rho, t, and eta

        for tstep in range(tstep0,tstep1,tstep_step): #each time step
            tstepstr=str(tstep).zfill(5)
            fname= '_'+field+tstepstr
            print(fname+'/'+str(tstep1),end='\r')

            # get binary data from file. change dtype if your output is a different accuracy
            data = np.fromfile(dir+fname, dtype=np.float64)
            data_end = data[:len(data)-y_dim*z_dim]
            data_begin = data[len(data)-y_dim*z_dim:]

            fig, ax = plt.subplots(2)
            fig.set_size_inches(16,7)
            ax[0].plot(range(len(data_end)), data_end)
            ax[1].plot(range(len(data_begin)), data_begin)
            if showPlot_bool:
                plt.show()
            #else:
            #    print('not showing this plot')
            plt.clf()
            plt.close()

            # read the data off the individual cores and put them into a single table
            data_cores = []
            for i in range(nz_cores*ny_cores):
                data_cores.append(data_begin[i*nz_dim*ny_dim:(i+1)*nz_dim*ny_dim])
            data_cores=np.array(data_cores)
            data_all=np.zeros([z_dim,y_dim])
            for i in range(nz_cores):
                for j in range(ny_cores):
                    cn = i*ny_cores+j
                    core_data=data_cores[cn]
                    for ii in range(nz_dim):
                        for jj in range(ny_dim):
                            da_idx_z = i*nz_dim +ii
                            da_idx_y = j*ny_dim +jj
                            data_all[da_idx_z][da_idx_y] = core_data[ii*ny_dim+jj]
                    
            #OPEN MESH FILE	
            Tablemesh=[]
            with open(meshdir) as f:
                reader = csv.reader(f)
                for col in reader:
                    Tablemesh.append(col)
                            

            y=[float(i) for i in Tablemesh[0][0:y_dim]] # get the radii coords
            z=[float(i) for i in Tablemesh[0][y_dim::2*y_dim]] # get the z coords
            Y=len(y)
            Z=len(z)
            arrY=[i for i in range(1,Y+1)]
            arrZ=[i for i in range(1,Z+1)]

            fig, ax = plt.subplots(2)
            fig.set_size_inches(16,7)
            [ax[0].scatter(y,data_all[i]) for i in range(0,z_dim,40)]
            [ax[1].scatter(data_all.T[i],z) for i in range(0,y_dim,10)]
            if showPlot_bool:
                plt.show()
            plt.clf()
            plt.close()

            horizontal_profiles=data_all
            vertical_profiles = data_all.T

            if writetocsv_bool:
                with open(dir+field+'_'+tstepstr+'.csv', 'w', newline='') as csvfile:
                    fwriter = csv.writer(csvfile, delimiter=' ')
                    for k in range(len(horizontal_profiles)):
                        fwriter.writerow(data_all[k])
                csvfile.close()
