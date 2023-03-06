# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 17:04:57 2021

@author: legraine
"""

import numpy as np 
import csv, math
from pylab import *
import matplotlib.pyplot as plt

#######################################################
#BestParamPlot chosen by the user (could be pasted directly from the ABR_kO_simulation.py program output)

#Best fit BestParamPlot
BestParamPlot = (0.03476178959337191, -0.05237833344209468, 0.19114578135745458, 0.7168699559694265, -0.004863292275220595, 7.630602781172232, 23.593183917258482, -0.33579051041153607, 12.19310038625117, 90.27137001778543, 6.430325029635019, -3.1122953720236968, 1232.519432640296, 21.25809026774471)


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

#######################################################
#Compute the derivative for each time step

def Phi(i,vt,parameters):
    
    if i/2>(2000-parameters[12]):
        test_threshold_gd = parameters[6]*Esi[i]+parameters[7]*Eco[i]+parameters[8]*EnOtr[i]+vt
        test_threshold_dg = parameters[6]*Esi[i]+parameters[7]*Eco[i]+parameters[8]*EnOtr[i]
    else :
        test_threshold_gd = parameters[6]*Esi[i]+parameters[7]*Eco[i]+parameters[13]*EnOtr[i]+vt
        test_threshold_dg = parameters[6]*Esi[i]+parameters[7]*Eco[i]+parameters[13]*EnOtr[i]
        
    if S[0]=="g":
        if test_threshold_gd>parameters[9] and test_threshold_dg>parameters[10]:
            S[0] = "d"
            S[1] = i/2*0.1
    else :
        if test_threshold_dg<parameters[10] and test_threshold_gd<parameters[9]:
            S[0] = "g"

    if S[0]=="g" :
        dvdt = -parameters[0]*Esitr[i]-parameters[1]*Ecotr[i]-parameters[2]*EnOtr[i]+parameters[3]
    else :
        dvdt = -parameters[0]*Esitr[i]-parameters[1]*Ecotr[i]-parameters[2]*EnOtr[i]+parameters[4]-vt/parameters[5]
            
    return(dvdt)

##########################################################
#Compute the modelled volume for the best BestParamPlot using the Runge–Kutta 4th order method

def modelledVolume(a,b,vi,n,BestParamPlot) :
    t = np.linspace(a,b,n+1)
    v.append(vi)
    state.append(S[1])
    pas = (b-a)/float(n)
    j = 0
    for k in range(n) :
        if S[0]=="g":
            state.append(0)
        else:
            state.append(1)
        k1 = Phi(2*k,v[k],BestParamPlot)
        k2 = Phi(2*k+1,v[k]+k1*pas/2.,BestParamPlot)
        k3 = Phi(2*k+1,v[k]+k2*pas/2.,BestParamPlot)
        k4 = Phi(2*k+2,v[k]+pas*k3,BestParamPlot)
        v.append(v[k]+pas/6.*(k1+2*k2+2*k3+k4))
    return v

##########################################################
#Best BestParamPlot

aEsi = BestParamPlot[0]
aEco = BestParamPlot[1]
aO = BestParamPlot[2]
ag = BestParamPlot[3]
ad = BestParamPlot[4]
tau = BestParamPlot[5]
kEsi = BestParamPlot[6]
kEco = BestParamPlot[7]
kO = BestParamPlot[8]
v0 = BestParamPlot[9]
v1 = BestParamPlot[10]
vi = BestParamPlot[11]
ageMPT = BestParamPlot[12]
agMPT = BestParamPlot[13]
S = ["g",0]


#Lists initialization
time=[];esinomega=[]; ecosomega=[] ; epsilon=[] ; O=[]
sea=[]
state=[]
Esi=[]; Eco=[];Esitr=[];Ecotr=[]
v=[];dvdt=[]
EnOtr=[]
EnO=[]

#####################################################################
#Recovery and processing of orbital data for the last 2Ma at the time step of 1ka

file = open('Orbital_parameters_2Ma_1ka.csv',"r")
filecsv = csv.reader(file, delimiter=';')

l = 0
for row in filecsv:
    time.append(float(row[0]))
    esinomega.append(float(row[1]))
    ecosomega.append(float(row[2]))
    O.append(float(row[3]))
    l=l+1
file.close()

#Normalization and truncation of BestParamPlot input
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

##########################################################
#Modelling of the ice volume for the best BestParamPlot fit

icevolume = modelledVolume(-2000,0,vi,2000,BestParamPlot)

#calcul de l'écart modele donnees a chaque pas de temps
residuals = []
sum_residuals = 0
for i in range (len(time)):
    sum_residuals = sum_residuals + (sea[i]-v[i])**2
    residuals.append(sea[i]-v[i])
    
    
##########################################################
#outputs data    
    
print("      Minimum residuals = " + str((sum_residuals)))   
print("      Average of residuals = "+ str(((((sum_residuals)/2000))**(1/2))))


##########################################################
#outputs figure
plt.close('all')
rcParams["figure.figsize"] = [20, 4]

#Comparison model-data for the best fit BestParamPlot
fig1, ax = plt.subplots()

ax.plot(time,sea,"r--", label="Data")
ax.plot(time,icevolume,'b',label="Model")
plt.gca().invert_yaxis()
xticks([0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
yticks([0,20,40,60,80,100,120])
xlim(2000,0)
ylim(-10,150)
plt.gca().invert_yaxis()
xlabel("Age (ka)",weight='bold')
ylabel("Ice volume (m sl)",weight='bold')

#state of the model for the best fit BestParamPlot
fig2, ax = plt.subplots()
  
ax.plot(time,state,"0.75", linewidth = 0.8 )
xticks([0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
xlim(0,2000)
yticks([0,1])
xlabel("Age (ka)",weight='bold')
ylabel("g or d (d=1, g=0)",weight='bold')
xlim(2000,0)

#residuals (model-data) for the best fit BestParamPlot
fig3, ax = plt.subplots()

ax.plot(time,residuals,'b--',label="Model")
show()
plt.gca().invert_yaxis()
xticks([0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
yticks([-40,-30,-20,-10,0,10,20,30,40])
xlim(2000,0)
xlabel("Age (ka)",weight='bold')
ylabel('Residuals (model-data)(msl)',weight='bold')
yticks([-40,-20,0,20,40])
