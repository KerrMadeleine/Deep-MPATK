import set_DeepMPATK_pars as par
import numpy as np;
import matplotlib.pyplot as plt;
import csv;

"""
This script reads a StagYY output file called rprof.dat to obtain the 
V_rms evolution data for use in determining "steady state"
TODO: remove the csv reading and just use np.loadtxt :P
"""

baseDir = par.STAGYY_OUTPUT_FOLDER
modeldirs = par.STAGYY_OUTPUT_MODS
Directories = [baseDir+modeldirs[i] for i in range(len(modeldirs))]
#Directories=Directories[par.STAGYY_OUTPUT_MODS_0:par.STAGYY_OUTPUT_MODS_1+1]
savefigure=par.READ_RPROF_SAVEFIG
showfigure=par.READ_RPROF_SHOWFIG

print(modeldirs)
print(Directories)
for i in range(par.STAGYY_OUTPUT_MODS_0,par.STAGYY_OUTPUT_MODS_1):
	print(i)
	Directory0 = Directories[i]+'/_rprof.dat'
	###

	###################################
	# Reference 
	###################################
	#0	z			#20	sprof2		#40 denprof2
	#1	Tmean		#21	sprof3 		#41	denprof3
	#2	Tmin		#22	wprof11		#42	airprof1
	#3	Tmax		#23	wprof21		#43	airprof2
	#4	vprof1		#24	wprof31		#44	airprof3
	#5	vprof2		#25	wprof12		#45	primprof1
	#6	vprof3		#26	wprof22		#46	primprof2
	#7	vzprof1		#27	wprof23		#47	primprof3
	#8	vzprof2		#28	dprof1		#48	ccprof1
	#9	vzprof3		#29	dprof2		#49	ccprof2
	#10	vhprof1		#30	dprof3		#50	ccprof3
	#11	vhprof2		#31	enprof1		#51	fprof1
	#12	vhprof3		#32	enprof2		#52	fprof2
	#13	etaprof1	#33	enprof3		#53	fprof3
	#14	etaprof2	#34	enprof4		#54	metalprof1
	#15	etaprof3	#35	enprof5		#55	metalprof2
	#16	eprof1		#36	cprof1		#56	metalprof3
	#17	eprof2		#37	cprof2		
	#18	eprof3		#38	cprof3
	#19	sprof1		#39	denprof1
	####################################


	# Don't go below this line unless you REally Wanna
	################################################

	TableData = []
	TimeSepData = []

	with open(Directory0) as f:
		reader = csv.reader(f, delimiter=" ", skipinitialspace=True)
		#headers = next(reader,None)
		#print(headers)
		for row in reader:
			TableData.append(row)
		
	TableData=TableData[1:] #get rid of blank first row

	temp=[]
	for row in TableData:
		
		if row[0][0]=='*':
			TimeSepData.append(temp)
			temp=[]
			temp.append(row[5])

		else:
			temp.append(row)

	TimeSepData=TimeSepData[1:]

	time=[]

	for el in TimeSepData:

		if len(el[0])==10:
			time.append(float(el[0])/31556952/10**6)

	Data0 = []
	temp2=[]
	for tstep in TimeSepData:
		data=tstep[1:]

		for ii in range(len(data)):
			temp2.append(data[ii])

		Data0.append(temp2)
		temp2=[]
	

	Data1=[]
	for j in range(len(time)):
		MomData=[]
		for k in range(57):

			parData=[]
			for row in Data0[j]:

				parData.append(float(row[k]))

			MomData.append(parData)

		Data1.append(MomData)
	
	rad=Data1[0][0] #initial depth list
	rad_km=[r/1000 for r in rad]
	
	ListAvgProfs=[]
	ViscProfs=[]
	VrmsProfs=[]
	ListAvgNames=[]
	ListTimes=[]
	Datakey=1
	Hrng=0
	last=1
	
	for ii in range(len(time))[0:-last:Hrng+1]:

		avgTprof=[np.mean([Data1[k][Datakey][j] for k in range(ii-Hrng,ii+Hrng+1)]) for j in range(0,len(rad))]
		viscprof=[np.mean([Data1[k][13][j] for k in range(ii-Hrng,ii+Hrng+1)]) for j in range(0,len(rad))]
		vrmsprof=[np.mean([Data1[k][4][j] for k in range(ii-Hrng,ii+Hrng+1)]) for j in range(0,len(rad))]
		
		ListAvgProfs.append(avgTprof)
		ViscProfs.append(viscprof)
		VrmsProfs.append(vrmsprof)
		viscprof=[]
		avgTprof=[]

		ListAvgNames.append(str(time[ii])+' Myrs')
		ListTimes.append(time[ii])
	
	L= len(ListAvgProfs)

	def mantle_vrms(r_km,vrmsprofiles,limr0,limr):

		dr_km=[(r_km[jj+1]-r_km[jj]) for jj in range(limr0,limr+1,1)]
		r_km=r_km[1:]
		vrms_mantle=[]

		for t in range(len(vrmsprofiles)):
			sum=0
			sum_ws=0

			for sh in range(len(dr_km)):
				sum=sum+dr_km[sh]*r_km[sh]*np.pi*(vrmsprofiles[t][limr0+sh])**2
				sum_ws=sum_ws+dr_km[sh]*r_km[sh]*np.pi

			vrms_mantle.append(np.sqrt(sum/sum_ws))

		return vrms_mantle
	##################################################

	oneeighth=len(VrmsProfs[0])//8
	r0=oneeighth*2
	r1=oneeighth*5
	vrms_0 = [VrmsProfs[ii][0] for ii in range(len(VrmsProfs))]
	vrms_1 = [VrmsProfs[ii][-1] for ii in range(len(VrmsProfs))]
	vmrs_mean=[np.mean(VrmsProfs[ii][0:r1]) for ii in range(len(VrmsProfs))]
	vrms_mean_stag=[np.mean(VrmsProfs[ii]) for ii in range(len(VrmsProfs))]
	rad_km.insert(0,0)

	weighted_vrms_mantle=mantle_vrms(rad_km,VrmsProfs,r0,r1)

	fig, ax = plt.subplots(1)
	fig.set_size_inches(13,5)
	ax.plot(vrms_0,color='blue',linewidth=3)
	ax.plot(vrms_1,color='orange',linewidth=3)
	ax.plot(vmrs_mean,linewidth=5,color='k')
	ax.plot(vrms_mean_stag,linewidth=5,color='red')
	ax.plot(weighted_vrms_mantle,linewidth=5,color='green')
	ax.grid()
	if savefigure:
		print(len(Directories))
		print(i)
		print(Directories)
		plt.savefig(Directories[i]+'a_VRMS_plot.png')
	#if showfigure:
	plt.show()
	plt.clf()
	plt.close()


	onset_file = open(par.CONV_ONSET_FILE+modeldirs[i][:-1]+'.txt', 'a')
	onset_file.write('\n Vrms:, '+ str(np.mean(weighted_vrms_mantle[30:])))
	onset_file.close()