#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  8 20:42:37 2020

@author: mallein

Comparison of detection methods for a closed facility.
"""

import matplotlib.pyplot as plt
import numpy.random as npr
from math import sqrt,exp,log, ceil, floor, erf

def phi(x,mu,sigma):
    """
    Cumulative distribution function for the normal distribution with mean mu and variance sigma^2

    Parameters
    ----------
    x : float
        Point value of the CDF.
    mu : float
        Mean of the normal distribution.
    sigma : float
        Standard deviation of the normal distribution.

    Returns
    -------
    float
        Cumulative distribution function of the normal distribution at point x.

    """
    return (1.0 + erf((x-mu) / (sigma*sqrt(2.0)))) / 2.0

def probDetection(mu,sigma,nb):
    """
    Probability of detecting one contaminated individual in a group of nb

    Parameters
    ----------
    mu : float
        Mean of the logarithm of the viral concentraiton.
    sigma : float
        Standard deviation of the logarithm of the viral concentration.
    nb : int
        Number of combined individuals in the sample. Assumed that nb-1 are non-contaminated.

    Returns
    -------
    float
        Probability of detecting the positive sample via standard rt-qPRC.

    """
    return phi(40-log(nb)/log(2),mu, sigma)

def pSampling(b,k,n):
    """
    Probability of choosing a contaminated individual in a group.

    Parameters
    ----------
    b : int
        Number of contaminated individuals.
    k : int
        Number of selected individuals.
    n : int
        Total number of individuals in a group.

    Returns
    -------
    float
        Probability that choosing k individuals among n, one of the b contaminated is sampled.

    """
    if k >=n-b:
        return 1
    else:
        eps = 1
        for j in range(k):
            eps = eps*(n-b-j)/(n-j)
        return 1-eps

def detectionImproved(nbInfectes,nbTests,nbSamplePerTest,nbIndividus,mu,sigma):
    """
    Simulation of the result of a PCR made on a random sample of a population with the following parameters.

    Parameters
    ----------
    nbInfectes : int
        Number of infected individuals in the population.
    nbTests : int
        Number of tests made in the population at the same time.
    nbSamplePerTest : int
        Number of individuals pooled together for a test.
    nbIndividus : int
        Total size of the population.
    mu : float
        Mean of the log of the Ct value
    sigma : float
        Std of the log of the Ct value

    Returns
    -------
    bool
        True if detected, False if undetected.

    """
    if nbTests > 0:
        if nbTests*nbSamplePerTest > nbIndividus:
            nbSamplePerTest = int(floor(nbIndividus/nbTests))
        listeConcentrations = [0 for j in range(nbTests)]
        listeIndividus = [nbSamplePerTest for j in range(nbTests)]
        for individual in range(nbInfectes):
            noIndiv = int(floor(npr.rand()*nbIndividus-individual))
            testNo=-1
            while noIndiv>=0 and testNo < nbTests-1:
                testNo+=1
                noIndiv-= listeIndividus[testNo]
            if noIndiv<0:
                listeIndividus[testNo] -=1
                listeConcentrations[testNo]+= 1
        detection = False
        for test in range(nbTests):
            detection = detection or listeConcentrations[test]/nbSamplePerTest> 0
        return detection
    else:
        return False

def runEpidemicsBis(timeWindow,rateInfection,asymptomaticRatio,totalNumber,timeSymp,nbTests,nbSamplePerTest,mu,sigma):
    """
    Run the onset of an epidemic in a community practicing regular screening tests.
    The epidemic runs until its first detection, whether thanks to a screening test
    returning positive, or an individual finally presenting symptoms.


    Parameters
    ----------
    timeWindow : float
        Time in between two screening tests. The epidemics starts at a time
        uniformly distributed between two tests.
    rateInfection : float
        Rate at which contamination occurs within the community.
    asymptomaticRatio : float
        Proportion of asymptomatic carriers within the community.
    totalNumber : int
        Number of individuals in this community.
    timeSymp : float
        Time before first detection of symptoms. DETERMINISTIC, BUT COULD BE MADE RANDOM.
    nbTests : int
        Number of tests used at each screening event. EACH TEST IS MADE INDEPENDENTLY FOR SIMPLIFICATION, SHOULD ADD SOME SYNERGY.
    nbSamplePerTest : int
        Number of individuals sampled when a test is made.
    probMesure : float
        Probability that the test will detect contamination if at least one contaminated individual belongs to it.

    Returns
    -------
    numberInfected : int
        Number of infected individuals when detection occurs.
    currentTime : float
        Time after infection at which detection occurs.

    """
    currentTime = 0
    numberInfected = 1
    nextDetection = npr.rand()*timeWindow
    firstSymptoms = 100000000000000000000000000000
    if npr.rand()>asymptomaticRatio:
        firstSymptoms = min(firstSymptoms, timeSymp+currentTime)
    
    isDetected = False
    while not isDetected:
        nextInfection = currentTime - log(npr.rand())/(rateInfection*numberInfected)
        nextEvent = min(nextInfection,nextDetection,firstSymptoms)
        if firstSymptoms == nextEvent:
            currentTime = firstSymptoms
            isDetected = True
        elif nextDetection == nextEvent:
            currentTime = nextDetection
            nextDetection+=timeWindow
            isDetected = detectionImproved(numberInfected, nbTests, nbSamplePerTest, totalNumber, mu, sigma)
        else : 
            currentTime = nextInfection
            numberInfected+=1
            if npr.rand() > asymptomaticRatio:
                firstSymptoms = min(firstSymptoms,timeSymp+currentTime)
    return numberInfected, currentTime

def detection(b,k,n):
    """
    Simulation of a perfect test in a population of n with b infected and k sampled.

    Parameters
    ----------
    b : int
        Number of infected.
    k : int
        Number of selected for the test.
    n : int
        Total size of the population.

    Returns
    -------
    bool
        True if detected (i.e. if one of the k selected is infected), False otherwise.

    """
    eps = 1
    for j in range(k):
        eps = eps*(n-b-j)/(n-j)
    return npr.rand() > eps

def runEpidemics(timeWindow,rateInfection,asymptomaticRatio,totalNumber,timeSymp,numberTested):
    """
    Run the onset of an epidemic in a community practicing regular perfect screening tests.
    The epidemic runs until its first detection, whether thanks to a screening test
    returning positive, or an individual finally presenting symptoms.

    Parameters
    ----------
    timeWindow : float
        Time in between two screening tests. The epidemics starts at a time
        uniformly distributed between two tests.
    rateInfection : float
        Rate at which contamination occurs within the community.
    asymptomaticRatio : float
        Proportion of asymptomatic carriers within the community.
    totalNumber : int
        Number of individuals in this community.
    timeSymp : float
        Time before first detection of symptoms.
    numberTested : int
        Number of individuals sampled when a test is made.

    Returns
    -------
    numberInfected : int
        Number of infected individuals when detection occurs.
    currentTime : float
        Time after infection at which detection occurs.

    """

    nextDetection = npr.rand()*timeWindow
    
    currentTime = 0
    numberInfected = 1
    firstSymptoms = totalNumber * timeSymp
    if npr.rand()>asymptomaticRatio:
        firstSymptoms = timeSymp+currentTime
    
    isDetected = False
    while not isDetected:
        nextInfection = currentTime-log(npr.rand())/(rateInfection*numberInfected)
        nextEvent = min(nextInfection,nextDetection,firstSymptoms)
        if nextEvent == nextDetection:
            currentTime = nextDetection
            nextDetection += timeWindow
            isDetected = detection(numberInfected,numberTested,totalNumber)
        elif nextEvent == firstSymptoms:
            isDetected = True
            currentTime = firstSymptoms
        else:
            currentTime = nextInfection
            numberInfected +=1
            if npr.rand() > asymptomaticRatio:
                firstSymptoms = min(firstSymptoms,currentTime+timeSymp)
    return numberInfected,currentTime

window = 2
rateInfection = 0.5
asymptomaticRatio = 0.3
totalNumber = 1000
timeSymp = 5
nbTry = 10000
sampleMax = 200
mu = 35
sigma = 4

res = []
resBis = []
for j in range(sampleMax):
    res.append(0)
    resBis.append(0)
    for trial in range(nbTry):
        a,b = runEpidemicsBis(window,rateInfection,asymptomaticRatio,totalNumber,timeSymp,2,j,mu,sigma)
        res[j]+=a
        resBis[j]+=b
    res[j] = res[j]/nbTry
    resBis[j] = resBis[j]/nbTry

plt.plot(res)
plt.show()
plt.clf()

plt.plot(resBis)
plt.show()
plt.clf()
