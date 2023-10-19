import numpy as np;
import matplotlib.pyplot as plt;
import csv;
import sys;
import ALLDATA4 as AD4

writeData=False

baseDir = '/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/fromPerl/'
modeldirs = ['1300_1/','1400_1/','1500_1/','1600_1/','1700_1/','1800_1/','1900_1/']
saveDir= '/Users/mkerr/VenusResearch/2023/STAGYY/Analysis/PLOTS/read_time/'
#baseDir = '/Users/mkerr/VenusResearch/2022/FinalModels/'
#modeldirs = ['t0_1600_100/','t0_1600_200/','t0_1600_300/','t0_1700_100/','t0_1700_200/','t0_1700_300/','t0_1800_100/', 't0_1800_200/', 't0_1800_300/' ]
Directories = [baseDir+modeldirs[i] for i in range(len(modeldirs))]
print(Directories)


RAEFF=AD4.ALLRA
Plottable=[]
Times=[]
FLOW=[]
FLOW_T=[]
FLOW_B=[]
VRMS = []
VMAX=[]
RA_eff=[]
NU_BOT=[]
NU_TOP=[]
for i in range(len(Directories)):
	print('the', i,'th', 'directory is:', modeldirs[i][:-1])
	mod_name = modeldirs[i][:-1]
	time_data_file = Directories[i]+'_time.dat'

	TableData = []
	istep= []
	time_yr = []
	F_top= []
	F_bot= []
	F_del =[]
	ra_eff = []
	nu_bot=[]
	nu_top=[]
	V_rms = []
	V_max=[]

	#with open('/Users/mkerr/VenusResearch/2022/FinalModels/15e25_17k3/_time.dat') as f:
	with open(time_data_file) as f:
		reader = csv.reader(f, delimiter=" ", skipinitialspace=True)

		for row in reader:
			TableData.append(row)

	Headers = TableData[0]
	print(Headers)
	noHeaderTable = TableData[1::] 
	Data=np.zeros((len(noHeaderTable),len(noHeaderTable[0])))
	for i in range(len(noHeaderTable)):
		for j in range(len(noHeaderTable[0])-1):
			if noHeaderTable[i][j]!='1.340780792994+154':
				Data[i][j]=noHeaderTable[i][j]
	for row in Data:
		istep.append(row[0])
		time_yr.append(row[1]/31556952/10**6)
		F_top.append(row[2]*1e3) #mW/m^2
		F_bot.append(row[3]*1e3) #mW/m^2
		ra_eff.append(row[13])
		nu_bot.append(row[15])
		nu_top.append(row[14])
		F_del.append(row[2]*1e3-row[3]*1e3) #mW/m^2
		V_rms.append(row[8])
		V_max.append(row[9])

	FLOW.append(F_del)
	FLOW_T.append(F_top)
	FLOW_B.append(F_bot)
	Times.append(time_yr)
	VRMS.append(V_rms)
	VMAX.append(V_max)
	RA_eff.append(ra_eff)
	NU_BOT.append(nu_bot)
	NU_TOP.append(nu_top)

	data = list(zip(time_yr,F_del))

	if writeData:
		np.savetxt('HF_mW-m2_'+mod_name+'.csv',data, delimiter=",")	
	
for i in range(len(Directories)):
	# FLOW
	fig, ax = plt.subplots(1)
	ax.plot(Times[i],FLOW[i], label=modeldirs[i][:-1],linewidth=3,color='green')
	ax.plot(Times[i],FLOW_T[i], label=modeldirs[i][:-1],linewidth=3,color='red')
	ax.plot(Times[i],FLOW_B[i], label=modeldirs[i][:-1],linewidth=3,color='blue')
	ax.set_xlim([0,Times[i][-1]+10])
	ax.set_ylim([-150,50])
	ax.set_xlabel('Time (Myrs)')
	ax.set_ylabel('$mW/m^2$')
	ax.set_title('$\Delta F$ over Time: '+modeldirs[i][:-1])
	ax.grid()
	#plt.savefig(saveDir+modeldirs[i][:-1]+'/HeatFlow.png')
	#plt.show()
	plt.clf()
	plt.close()
	
	#VRMS
	fig, ax = plt.subplots(1)
	print(VRMS[i][-20:])
	ax.plot(Times[i],VRMS[i], label=modeldirs[i][:-1],linewidth=3)
	ax.plot(Times[i],VMAX[i], label=modeldirs[i][:-1],linewidth=3)
	ax.set_xlim([0,Times[i][-1]+10])
	ax.set_ylim([0,1e-8])
	ax.set_xlabel('Time (Myrs)')
	ax.set_ylabel('m/s')
	ax.set_title('$V_{rms}$ over Time: '+modeldirs[i][:-1])
	ax.grid()
	#plt.savefig(saveDir+modeldirs[i][:-1]+'/Vrms.png')
	#plt.show()
	plt.clf()
	plt.close()

	# RA_EFF
	fig, ax = plt.subplots(1)
	ax.plot(Times[i],RA_eff[i], label=modeldirs[i][:-1],linewidth=3)
	ax.set_xlim([0,Times[i][-1]+10])
	#ax.set_ylim([0,200])
	ax.set_xlabel('Time (Myrs)')
	ax.set_ylabel('Ra_eff')
	ax.set_title('RA_EFF over Time: '+ modeldirs[i][:-1])
	ax.grid()
	#plt.savefig(saveDir+modeldirs[i][:-1]+'/ra_eff.png')
	#plt.show()
	plt.clf()
	plt.close()

	# NU BOT
	fig, ax = plt.subplots(1)
	ax.plot(Times[i],NU_BOT[i], label=modeldirs[i][:-1],linewidth=3)
	ax.set_xlim([0,Times[i][-1]+10])
	#ax.set_ylim([0,200])
	ax.set_xlabel('Time (Myrs)')
	ax.set_ylabel('Nu bot')
	ax.set_title('NU_BOT: '+ modeldirs[i][:-1])
	ax.grid()
	#plt.savefig(saveDir+modeldirs[i][:-1]+'/nu_bot.png')
	#plt.show()
	plt.clf()
	plt.close()

	# NU TOP
	fig, ax = plt.subplots(1)
	ax.plot(Times[i],NU_TOP[i], label=modeldirs[i][:-1],linewidth=3)
	ax.set_xlim([0,Times[i][-1]+10])
	#ax.set_ylim([0,200])
	ax.set_xlabel('Time (Myrs)')
	ax.set_ylabel('Nu top')
	ax.set_title('NU_TOP: '+ modeldirs[i][:-1])
	ax.grid()
	#plt.savefig(saveDir+modeldirs[i][:-1]+'/nu_top.png')
	#plt.show()
	plt.clf()
	plt.close()






fig, ax = plt.subplots(1)
for i in range(len(Directories)):
	ts=len(RAEFF[i])
	ax.plot(RAEFF[i],NU_BOT[i][0:40*ts:40], label=modeldirs[i][:-1],linewidth=1)
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid(which='major', color='k', linestyle='-')
ax.grid(which='minor', color='gray', linestyle='--')
#plt.savefig(saveDir+modeldirs[i][:-1]+'/nu_top.png')
plt.show()
plt.clf()
plt.close()

###############################################################
#### REFERENCE FOR TIME FILE #########
###############################################################
# 0 : istep (time step)
# 1 : time  (time in seconds)
# 2 : F_top (heat flow from the top)
# 3 : F_bot (heat flow from the bottom)
# 4 : Tmin  (min Temp)
# 5 : Tmean (avg Temp)
# 6 : Tmax (max Temp)
# 7 : Vmin 
# 8 : Vrms 
# 9 : Vmax  
# 10: eta_min 
# 11: eta_mean 
# 12: eta_max 
# 13: ra_eff  
# 14: Nu_top
# 15: Nu_bot
# 16: C_min
# 17: C_mean
# 18: C_max
# 19: F_mean
# 20: F_max
# 21: erupt_rate
# 22: erupta
# 23: erupt_heatflux
# 24:entrainment
# 25: Cmass_error
# 26: H_int
###############################################################


