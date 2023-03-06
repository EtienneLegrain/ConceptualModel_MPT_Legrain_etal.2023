# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 12:34:42 2021

@author: legraine
"""

import numpy as np 
import csv, math
import emcee
from multiprocessing import Pool
from pylab import *

#####################################################################
#Parameters chosen by the user

#Define the approximate initial position from which is computed each start position of the n walkers
StartPosition = (0.03, -0.2, 1, 0.02, 14, 2, 2, -7, 11, 80, -4, 10)




#Number of walkers (verifying ; nwalkers > 2 * number of parameters)
nwalkers = 30

#Number of iterations 
nbiterations = 10000

#Define the first position of each walkers relatively to StartPosition. When walkers_jump is high, walkers are far from StartPosition.
walkers_jump = 0.1

#####################################################################
#Compute the mean

def mean(sample) :
    size = len(sample)
    mean = sum(sample)/size
    return mean

#####################################################################
#Compute the standard deviation

def stdev(sample) :
    n = len(sample) 
    mq = mean(sample)**2
    s = sum([x**2 for x in sample])
    variance = s/n-mq
    stdeviation = math.sqrt(variance)
    return stdeviation

#####################################################################
#Normalise a distribution

def normalise(sample):
    norm = []
    average = mean(sample)
    stdeviation = stdev(sample)
    for i in range(len(sample)):
        norm.append((sample[i]-average)/stdeviation)
    return norm

#####################################################################
#Apply the truncation function from Parrenin and Paillard, 2012 

def truncation(sample):
    trunc = []
    for i in range(len(sample)):
        if sample[i]<=0:
            trunc.append(sample[i]+math.sqrt(4*pow(1.06587,2)+pow(sample[i],2))-2*1.06587)
        else:
            trunc.append(sample[i])
    return trunc

#####################################################################
#Interpolate to artificially increase resolution (multiply by fac the number of points)

def interpol(sample, fac):
    new_sample = []
    tab = range(fac)
    for j in range(len(sample)):
        difference = sample[j+1]-sample[j]
        new_difference = difference/fac
        for x in tab:
            new_sample.append(sample[j]+x*new_difference)
        if j == len(sample)-2:
            break
    new_sample.append(sample[len(sample)-1]) #specific case for the last point
    return new_sample

#####################################################################
#Compute the derivative for each time step

def Phi(i,vt,parameters):
    test_threshold_gd = parameters[6]*Esi[i]+parameters[7]*Eco[i]+parameters[8]*EnOtr[i]+vt
    test_threshold_dg = parameters[6]*Esi[i]+parameters[7]*Eco[i]+parameters[8]*EnOtr[i]
    if S[0] == "g":
        if test_threshold_gd > parameters[9] and test_threshold_dg > parameters[10]:
            S[0] = "d"
            S[1] = i/2*0.1
    else :
        if test_threshold_dg < parameters[10] and test_threshold_gd < parameters[9]:
            S[0] = "g"
    if S[0] == "g" :
        dvdt = -parameters[0]*Esitr[i]-parameters[1]*Ecotr[i]-parameters[2]*EnOtr[i]+parameters[3]
    else :
        dvdt = -parameters[0]*Esitr[i]-parameters[1]*Ecotr[i]-parameters[2]*EnOtr[i]+parameters[4]-vt/parameters[5]
    return dvdt

#####################################################################
#Compute the modelled volume for a set of input parameters using the Runge–Kutta 4th order method

def modelledVolume(a,b,vi,n,parameters):
    vt = np.zeros(2001)
    vt[0] = vi
    step = (b-a)/float(n)
    for k in range(n) :
        k1 = Phi(2*k,vt[k],parameters)
        k2 = Phi(2*k+1,vt[k]+k1*step/2.,parameters)
        k3 = Phi(2*k+1,vt[k]+k2*step/2.,parameters)
        k4 = Phi(2*k+2,vt[k]+step*k3,parameters)
        vt[k+1]=vt[k]+step/6.*(k1+2*k2+2*k3+k4)
    return vt   

#####################################################################
#Compute the residuals between model and data for a set of input parameters 

def cost_function_negative(parameters):
    residuals = 0
    Vtotal = modelledVolume(-2000,0,parameters[11],2000,parameters)
    delta = (Vtotal-sea)**2
    residuals = -np.sum(delta)
    return residuals

#####################################################################
#Lists initialization

time=[]; esinomega=[]; ecosomega=[] ; epsilon=[] ; O=[]; Esi=[]; Eco=[]; Esitr=[]; Ecotr=[]; sea=np.array([])

#####################################################################
#Recovery and processing of orbital data for the last 2Ma at the time step of 1ka

file = open('Orbital_Parameters_2Ma_1ka_synthetic.csv',"r")
filecsv = csv.reader(file, delimiter=';')

l = 0
for row in filecsv:
    time.append(float(row[0]))
    esinomega.append(float(row[1]))
    ecosomega.append(float(row[2]))
    O.append(float(row[3]))
    l=l+1
file.close()

#Normalization and truncation of parameters input
EnO = normalise(O)
Esi = normalise(esinomega)
Eco = normalise(ecosomega)
Esitr = truncation(Esi)
Ecotr = truncation(Eco)
EnOtr = truncation(EnO)
Esitr = normalise(Esitr)
Ecotr = normalise(Ecotr)
EnOtr = normalise(EnOtr)

#Interpolation to get data at the time step of 500 years (for half-step Runge-Kutta computation)
Esi = interpol(Esi,2)
Eco = interpol(Eco,2)
Esitr = interpol(Esitr,2)
Ecotr = interpol(Ecotr,2)
EnOtr = interpol(EnOtr,2)

#####################################################################
#Recovery and processing of ice volume data from Berends et al., 2021 for the last 2Ma at the time step of 1ka

file = open('Sea-level_Berends_2Ma_1ka.csv',"r")#données
filecsv = csv.reader(file, delimiter=';')

k = 0
for row in filecsv:
    if k>2000:
        k=k+1
    else:
        sea=np.append(sea, float(row[0]))
        k=k+1
file.close()

#####################################################################
#Assign to each model parameter the corresping value of the StartPosition input 

aEsii = StartPosition[0]
aEcoi = StartPosition[1]
aOi = StartPosition[2]
agi = StartPosition[3]
adi = StartPosition[4]
taui = StartPosition[5]
kEsii = StartPosition[6]
kEcoi = StartPosition[7]
kOi = StartPosition[8]
v0i = StartPosition[9]
v1i = StartPosition[10]
vii = StartPosition[11]
S = ["g",0]

param = [aEsii,aEcoi,aOi,agi,adi,taui,kEsii,kEcoi,kOi,v0i,v1i,vii]


#####################################################################
#Exploration of the parameters space using a markov chain methodo (MCMC) coupled with a random walk at n walkers using the eemc hammer (Foreman-Mackey, 2013)


#Define the initial position of each walkers from StartPosition and walkers_jump input values
ndim = len(param)
WalkersIni = np.zeros((nwalkers, ndim))
StdDevParam = np.zeros(len(param))
for i in range (len(param)):
    StdDevParam[i] = param[i]*walkers_jump
        
for j in range (nwalkers):
    for i in range (len(param)):
        WalkersPos = param[i] + np.random.normal(0,abs(StdDevParam[i]))
        WalkersIni[j][i] = WalkersPos

################# using a Linux or Apple machine (parallel hearts computation)

#with Pool() as pool:
#One iteration computation 
#    steps = emcee.EnsembleSampler(nwalkers, ndim, cost_function_negative, pool = pool)
#Compute from the intial position of walkers WalkersIni n times with n correspond to nbiterations
#    steps.run_mcmc(WalkersIni,nbiterations)

################# using a Windows machine

#One iteration computation 
steps = emcee.EnsembleSampler(nwalkers, ndim, cost_function_negative)

#Compute from the intial position of walkers WalkersIni n times with n correspond to nbiterations
steps.run_mcmc(WalkersIni,nbiterations)

#####################################################################
#Extraction of the best parameters list to copy and paste in the SYN_simulation_plot.py program

#Recovering of parameters from which we obtain the minimum residuals
maxindex = np.argmax(steps.get_log_prob(flat=True))
meilleurproba = np.max(steps.get_log_prob(flat=True))
maxfit = steps.get_chain(flat=True)[maxindex]
variables = maxfit

#Convert variables into a list
BestParam = []
for i in range (len(variables)):
    BestParam.append(variables[i])
    
print("Minimum residuals = " + str(abs(meilleurproba)))   
print ("Average of residuals = "+ str((abs(int(meilleurproba)/2000))**(1/2)))
print ("Best fit parameters are : ", BestParam)    

