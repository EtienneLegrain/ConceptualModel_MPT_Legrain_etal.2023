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
#Parameters chosen by the user (could be pasted directly from the SYN_simulation.py program output)

#Best fit parameters
BestParamPlot = (0.030033212527457498, -0.18869940707356844, 1.3952897980508958, 0.022590169184870553, 14.19834754619912, 2.345975160662065, 1.7955452379565073, -7.489542573035756, 11.163290813156008, 80.43108662202502, -4.027936921282937, 10.697790966594265)


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

def Phi(i,v):
    test_treshold_gd = kEsi*Esi[i]+kEco*Eco[i]+kO*EnOtr[i]+v
    test_treshold_dg = kEsi*Esi[i]+kEco*Eco[i]+kO*EnOtr[i]
    if S[0]=="g":
        if test_treshold_gd > v0 and test_treshold_dg > v1:
            S[0] = "d"
            S[1] = i/2
    else :
        if test_treshold_dg < v1 and test_treshold_gd < v0:
            S[0] = "g"
            Term_duration = i/2-S[1] #Compute the duration of a termination
            Term_start = 2000-S[1] #Compute the start of a termination
            ListDuration.append(Term_duration)
            ListStart.append(Term_start)

    if S[0]=="g" :
        dvdt = -aEsi*Esitr[i]-aEco*Ecotr[i]-aO*EnOtr[i]+ag
    else :
        dvdt = -aEsi*Esitr[i]-aEco*Ecotr[i]-aO*EnOtr[i]+ad-v/tau
    return dvdt

##########################################################
#Compute the modelled volume for the best parameters using the Runge–Kutta 4th order method

def modelledVolume(a,b,vi,n) :
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
        k1 = Phi(2*k,v[k])
        k2 = Phi(2*k+1,v[k]+k1*pas/2.)
        k3 = Phi(2*k+1,v[k]+k2*pas/2.)
        k4 = Phi(2*k+2,v[k]+pas*k3)
        v.append(v[k]+pas/6.*(k1+2*k2+2*k3+k4))
    return v

##########################################################
#Best parameters

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
S = ["g",0]


#Lists initialization
time=[];esinomega=[]; ecosomega=[] ; epsilon=[] ; O=[]
sea=[]
state=[]
Esi=[]; Eco=[];Esitr=[];Ecotr=[]
v=[];dvdt=[]
EnOtr=[]
EnO=[]
ListDuration=[]
ListStart=[]

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

##########################################################
#Modelling of the ice volume for the best parameters fit

icevolume = modelledVolume(-2000,0,vi,2000)

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
print("      Termination duration  : " + str(ListDuration))
print("      Start of termination : " + str(ListStart))

##########################################################
#outputs figure
plt.close('all')
rcParams["figure.figsize"] = [20, 4]

#Comparison model-data for the best fit BestParamPlot
fig1, ax = plt.subplots()

ax.plot(time,sea,"r--", label="Data")
ax.plot(time,icevolume,'k',label="Model")
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

ax.plot(time,residuals,'k--',label="Model")
show()
plt.gca().invert_yaxis()
xticks([0,100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000])
yticks([-40,-30,-20,-10,0,10,20,30,40])
xlim(2000,0)
xlabel("Age (ka)",weight='bold')
ylabel('Residuals (model-data)(msl)',weight='bold')
yticks([-40,-20,0,20,40])