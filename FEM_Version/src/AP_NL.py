import readxtor as rx
import plotxtor as px
import numpy as np
from numpy.polynomial import polynomial
import matplotlib.pyplot as plt
import random

def calculate_AP_NL(psi, P):
    # Input parameters
    npsi = len(psi)
    
    # Output parameter
    AP_NL = np.zeros(10)
    
    # Local variables
    Pprime = np.zeros(npsi)
    
    # Calculate derivative of P with respect to psi
    for i in range(1, npsi-1):
        Pprime[i] = (P[i+1] - P[i-1]) / (psi[i+1] - psi[i-1])
        
    # Handle edge cases for Pprime at psi[0] and psi[npsi-1]
    Pprime[0] = (P[1] - P[0]) / (psi[1] - psi[0])
    Pprime[npsi-1] = (P[npsi-1] - P[npsi-2]) / (psi[npsi-1] - psi[npsi-2])
    
    # # Fit Pprime to a 9th degree polynomial using polyfit
    # AP_NL = polynomial.polyfit(psi, Pprime, deg=9)
    
    # Randomly split data into training and testing sets
    indices = list(range(npsi))
    random.shuffle(indices)
    split_index = int(0.9 * npsi)  # Split at 80% of the data
    train_indices, test_indices = indices[:split_index], indices[split_index:]
    psi_train, psi_test = psi[train_indices], psi[test_indices]
    Pprime_train, Pprime_test = Pprime[train_indices], Pprime[test_indices]
    
    # Fit Pprime to a 9th degree polynomial using polyfit on the training set
    AP_NL = polynomial.polyfit(psi_train, Pprime_train, deg=9)

    # Evaluate the fitted polynomial on the testing set
    Pprime_test_fit = polynomial.polyval(psi_test, AP_NL)
    fit_error = np.linalg.norm(Pprime_test-Pprime_test_fit)/np.linalg.norm(Pprime_test)
    
    # Plot Pprime and its fit
    plt.figure()
    plt.clf()
    plt.plot(psi,Pprime,lw=2, label='original')
    plt.plot(psi,np.polyvanp.polyval(AP_NL,psi),(AP_NL[::-1],psi[::-1]),lw=2, label='fit')
    plt.xlabel('Psi'); plt.ylabel('Pprime')
    plt.legend(loc='best')

    return AP_NL, fit_error

PsiM = rx.RawReadField('FRC','PsiMesh')
P = rx.RawReadField('FRC', 'Pressure')
J = rx.RawReadField('FRC','JacobMesh')
Vprime = np.mean(J[:,:-1],axis=1)

PV = P*Vprime**(5.0/3.0)
PV[-1] = 0

P_oldV_old = np.load('P_oldV_old.npy')
tol = 1e-4
error = np.linalg.norm(PV - P_oldV_old) / np.linalg.norm(P_oldV_old)
print(error)

if error > tol:
    P = P_oldV_old/Vprime**(5.0/3.0)
    P[-1] = 0
    AP_NL,fit_error = calculate_AP_NL(PsiM, P)
    print(AP_NL)
    print('Pprime fit error: {}'.format(fit_error))
else:
    print('error < tol')
