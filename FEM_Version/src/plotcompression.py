import readxtor as rx
import plotxtor as px

import numpy as np
import matplotlib.pyplot as plt

Z_060 = rx.RawReadField('Data/FRC_0.60','ZZ'); Z_070 = rx.RawReadField('Data/FRC_0.70','ZZ')
R_060 = rx.RawReadField('Data/FRC_0.60','RR'); R_070 = rx.RawReadField('Data/FRC_0.70','RR')
psi_060 = rx.RawReadField('Data/FRC_0.60','psi2d'); psi_070 = rx.RawReadField('Data/FRC_0.70','psi2d')

ZM_060 = rx.RawReadField('Data/FRC_0.60', 'ZMesh'); ZM_070 = rx.RawReadField('Data/FRC_0.70', 'ZMesh')
RM_060 = rx.RawReadField('Data/FRC_0.60', 'RMesh'); RM_070 = rx.RawReadField('Data/FRC_0.70', 'RMesh') 
PsiM_060 = rx.RawReadField('Data/FRC_0.60','PsiMesh'); PsiM_070 = rx.RawReadField('Data/FRC_0.70','PsiMesh')

P_060 = rx.RawReadField('Data/FRC_0.60', 'Pressure'); P_070 = rx.RawReadField('Data/FRC_0.70', 'Pressure')
S_060 = rx.RawReadField('Data/FRC_0.60','S'); S_070 = rx.RawReadField('Data/FRC_0.70','S')
J_060 = rx.RawReadField('Data/FRC_0.60','JacobMesh'); J_070 = rx.RawReadField('Data/FRC_0.70','JacobMesh')
Vprime_060 = np.mean(J_060[:,:-1],axis=1); Vprime_070 = np.mean(J_070[:,:-1],axis=1)

# px.Figure(Z_060,R_060,psi_060,equal=False,fignumber=1)
# plt.xlabel('Z'); plt.ylabel('R')

# px.Figure(Z_070,R_070,psi_070,equal=False,fignumber=2)
# plt.xlabel('Z'); plt.ylabel('R')

# npsi_060, ntheta_060 = ZM_060.shape 
# PsiM_060_2D = np.reshape(np.repeat(PsiM_060,ntheta_060),(npsi_060,ntheta_060)) 
# px.Figure(ZM_060,RM_060,PsiM_060_2D,equal=False,fignumber=3)
# plt.xlabel('Z'); plt.ylabel('R')

# npsi_070, ntheta_070 = ZM_070.shape 
# PsiM_070_2D = np.reshape(np.repeat(PsiM_070,ntheta_070),(npsi_070,ntheta_070))
# px.Figure(ZM_070,RM_070,PsiM_070_2D,equal=False,fignumber=4)
# plt.xlabel('Z'); plt.ylabel('R')

plt.figure(5); plt.clf()
px.Plot(S_060,P_060*Vprime_060**(5/3),lw=2,label='psiedge = -0.60')
px.Plot(S_070,P_070*Vprime_070**(5/3),lw=2,label='psiedge = -0.70')
plt.xlabel('$S$'); plt.ylabel('$PV\'^{5/3}$')
plt.legend()
