import numpy as np
import matplotlib.pyplot as plt

# ---------------- Plot viscosity --------------------------
T_model = np.linspace(1200, 2100, 100)
T_LB = np.linspace(18,80,100)

A = 7.72
B = -0.18
C = 0.00087

R_gas=8.314462618
E_act = 240e3

A_eta = 1e21/np.exp(E_act/(R_gas * 1600))

eta_LB = np.exp(A+B*T_LB+C*T_LB**2)

eta_model = A_eta * np.exp(E_act/(R_gas * T_model))

fig, ax = plt.subplots(1,2)
fig.set_size_inches(5,3)
ax[0].plot(T_LB, eta_LB/np.min(eta_LB), color='b')
ax[1].plot(T_model, eta_model/np.min(eta_model),color='k')
ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].grid()
ax[1].grid()
plt.show()
plt.clf()
plt.close()
