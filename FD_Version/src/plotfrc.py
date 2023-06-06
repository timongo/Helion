import readxtor as rx
import plotxtor as px
import numpy as np
import matplotlib.pyplot as plt

Z = rx.RawReadField('FRC','ZZ')
R = rx.RawReadField('FRC','RR')
psi = rx.RawReadField('FRC','psi')
plt.figure(1)
plt.clf()
px.Figure(Z,R,psi,equal=False)
plt.xlabel('Z'); plt.ylabel('R')

P = rx.RawReadField('FRC','Pressure')
S = rx.RawReadField('FRC','s')
J = rx.RawReadField('FRC','JacobMesh')
Vprime = np.mean(J[:,:-1],axis=1)
plt.figure(2)
plt.clf()
px.Plot(s,P*Vprime**(5/3),lw=2)
plt.xlabel('$S$'); plt.ylabel('$PV\'^{5/3}$')
