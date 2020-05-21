#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 09:00:07 2020

@author: mallein

Copy of TestBayesien for the purpose of testing
"""

import numpy as np
import numpy.random as npr
from math import log,exp,ceil, sqrt, floor
import matplotlib.pyplot as plt


def normalize(liste):
    """
    Normalize the list and return the median, 5% and 95% quantile
    
    Parameters
    ----------
    liste : float list
        .

    Returns
    -------
    med : float
        Index of the median.
    bs : int
        Index of the 95% quantile.
    bi : int
        Index of the 5% quantile.

    """
    med = 0
    bs = 0
    bi = 0
    s = sum(liste)
    tot = 0
    for k in range(len(liste)):
        liste[k] = liste[k]/s
        tot +=liste[k]
        if bi==0 and tot > .05:
            bi= k-1
        if med==0 and tot > .5:
            med= (2*k-1)/2
        if bs==0 and tot > .95:
            bs= k
    return med,bs,bi

def chooseSize(liste,step):
    """
    Choose the optimal size for the group testing depending on the current knowledge of the law

    Parameters
    ----------
    liste : float
        Law of the prior, summing to 1.
    stezp : float
        Step size of the discretization of the law of the prior.

    Returns
    -------
    N : int
        Optimal choice for the size of the group.

    """
    l = len(liste)
    thetaest = 0
    for k in range(l):
        thetaest+=step*liste[k]*k
    return int(ceil(log(0.2)/log(1 - thetaest)))



def measureOfPrevalenceUnidim(prevalence, nombreTests,Nmax,step):
    """
    Simulation of a Bayesian measure of prevalence in a population using Group testing

    Parameters
    ----------
    prevalence : float
        Probability for an individual to be contaminated.
    nombreTests : int
        Number of tests made during the experiment.
    Nmax : int
        Maximal number of individuals to sample in a single group
    step : float
        Width of the grid on which the prior is taken

    Returns
    -------
    listeMed : float list
        Median of the prior as a function of number of tests made
    listeBi : float list
        5% quantile of the prior as a function of number of tests made.
    listeBs : float list
        95% quantile of the prior as a function of number of tests made.

    """
    listetests = [j+1 for j in range(nombreTests)]
    j=0
    loiPrior = []
    while j*step < 1:
        theta = j*step
        loiPrior.append(theta*(1-theta))
        j+=1
    nbPoints = j
    med, bs, bi = normalize(loiPrior)
    
    listeMed = [med*step]
    listeBi = [bi*step]
    listeBs = [bs*step]
    
    for test in listetests:
        N = chooseSize(loiPrior,step)
        N = min(N,Nmax)

        resTest = 1*(npr.rand() > (1-prevalence)**N)
        for k in range(nbPoints):
            theta = k*step
            eps = (1-theta)**N
            if resTest == 1:
                loiPrior[k] = loiPrior[k]*(1-eps)
            else:
                loiPrior[k] = loiPrior[k]*eps
        med,bs,bi = normalize(loiPrior)
        listeMed.append(med*step)
        listeBi.append(bi*step)
        listeBs.append(bs*step)

    plt.plot(listeMed)
    plt.plot(listeBi)
    plt.plot(listeBs)
    plt.savefig("images/copymeasureOfPrevalenceUniDim.pdf")
    plt.show()
    plt.clf()
    fichier = open('csvForImage/measureOfPrevalenceUniDim.csv','w')
    fichier.write("Evolution of the estimation of the prevalence {} with maximal group size {}".format(prevalence,Nmax))
    fichier.write("Median : ")
    for el in listeMed:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Inf : ")
    for el in listeBi:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Sup : ")
    for el in listeBs:
        fichier.write(str(el)+" , ")
    fichier.write('\n')
    fichier.close()

    return listeMed, listeBi, listeBs

def binormalize(liste, ratio):
    """
    Normalize the 2-D list and return the median, 5% and 95% quantile of each direction as well as the combined version
    
    Parameters
    ----------
    liste : float list list
    ratio : float
        Proportion of individuals in the first coordinate.

    Returns
    -------
    med : float
        Index of the median.
    bs : int
        Index of the 95% quantile.
    bi : int
        Index of the 5% quantile.
    med1 : float
        Index of the median of the first subpopulation.
    bs1 : int
        Index of the 95% quantile of the first subpopulation.
    bi1 : int
        Index of the 5% quantile of the first subpopulation.
    med2 : float
        Index of the median of the second subpopulation.
    bs2 : int
        Index of the 95% quantile of the second subpopulation.
    bi2 : int
        Index of the 5% quantile of the second subpopulation.

    """
    s = 0
    l = len(liste)
    for j in range(l):
        s += sum(liste[j])
    listeAve = [0 for j in range(l)]
    liste1 = [0 for j in range(l)]
    liste2 = [0 for j in range(l)]
    
    for j in range(l):
        for k in range(l):
            liste[j][k] = liste[j][k]/s
            liste1[j] += liste[j][k]
            liste2[k] += liste[j][k]
            listeAve[int(floor(ratio*j+(1-ratio)*k))] += liste[j][k]
    
    med, bs, bi = normalize(listeAve)
    med1, bs1, bi1 = normalize(liste1)
    med2, bs2, bi2 = normalize(liste2)

    return med,bs,bi, med1, bs1, bi1, med2, bs2, bi2, listeAve


def measureOfPrevalenceBidim(prev1, prev2, ratio, nombreTests,Nmax,step):
    """
    Simulation of a Bayesian measure of prevalence in a population using Group testing.
    We assume a non-homogeneous population consisting of two subpopulations with different prevalences.

    Parameters
    ----------
    prev1 : float
        Probability for an individual to be contaminated in the first subpopulation.
    prev2 : float
        Probability for an individual to be contaminated in the second subpopulation.
    ratio : float
        Proportion of individuals belonging to the first subpopulation.
    nombreTests : int
        Number of tests made during the experiment.
    Nmax : int
        Maximal number of individuals to sample in a single group
    step : float
        Width of the grid on which the prior is taken

    Returns
    -------
    listeMed : float list
        Median of the prior as a function of number of tests made
    listeBi : float list
        5% quantile of the prior as a function of number of tests made.
    listeBs : float list
        95% quantile of the prior as a function of number of tests made.
    listeMed1 : float list
        Median of the prior as a function of number of tests made in the first subpopulation.
    listeBi1 : float list
        5% quantile of the prior as a function of number of tests made in the first subpopulation.
    listeBs1 : float list
        95% quantile of the prior as a function of number of tests made in the first subpopulation.
    listeMed2 : float list
        Median of the prior as a function of number of tests made in the second subpopulation.
    listeBi2 : float list
        5% quantile of the prior as a function of number of tests made in the second subpopulation.
    listeBs2 : float list
        95% quantile of the prior as a function of number of tests made in the second subpopulation.

    """
    listetests = [j+1 for j in range(nombreTests)]
    numberPoints=floor(1/step)+1
    loiPrior = [ [j*step*(1-j*step)*k*step*(1-k*step) for k in range(numberPoints)] for j in range(numberPoints)]
    
    med, bs, bi, med1, bs1, bi1, med2, bs2, bi2, loiAve = binormalize(loiPrior,ratio)
    
    listeMed = [med*step]
    listeBi = [bi*step]
    listeBs = [bs*step]
    listeMed1 = [med1*step]
    listeBi1 = [bi1*step]
    listeBs1 = [bs1*step]
    listeMed2 = [med2*step]
    listeBi2 = [bi2*step]
    listeBs2 = [bs2*step]
    
    for test in listetests:
        print(test)
        N = chooseSize(loiAve,step)
        N = min(N,Nmax)
        if test < 1000:
            N1= 0
            for k in range(N):
                if npr.rand()<ratio:
                    N1+=1
            N2 = N - N1
            probaTestNeg = ((1-prev1)**N1)* ((1-prev2)**N2)
            resTest = 1*(npr.rand() > probaTestNeg)
            for j in range(numberPoints):
                for k in range(numberPoints):
                    theta = j*step
                    phi = k*step
                    eps = ((1-theta)**N1)*((1-phi)**N2)
                    if resTest == 1:
                        loiPrior[j][k] = loiPrior[j][k]*(1-eps)
                    else:
                        loiPrior[j][k] = loiPrior[j][k]*eps
        else:
            print(test)
            loiAve1 = [sum(loiPrior[j]) for j in range(numberPoints)]
            loiAve2 = [sum(loiPrior[:][j]) for j in range(numberPoints)]
            N1= chooseSize(loiAve1,step)
            N2= chooseSize(loiAve2,step)
            probaTest1Neg = (1 - prev1)**N1
            probaTest2Neg = (1-prev2)**N2
            resTest1 = 1*(npr.rand()>probaTest1Neg)
            resTest2 = 1*(npr.rand()>probaTest2Neg)
            for j in range(numberPoints):
                for k in range(numberPoints):
                    theta = j * step
                    phi = k * step
                    eps1 = (1-theta)**N1
                    eps2 = (1-phi)**N2
                    if resTest1 == 1:
                        loiPrior[j][k] = loiPrior[j][k]*(1-eps1)
                    else:
                        loiPrior[j][k] = loiPrior[j][k]*eps1
                    if resTest2 == 1:
                        loiPrior[j][k] = loiPrior[j][k]*(1-eps2)
                    else:
                        loiPrior[j][k] = loiPrior[j][k]*eps2

        med, bs, bi, med1, bs1, bi1, med2, bs2, bi2, loiAve = binormalize(loiPrior,ratio)

        listeMed.append(med*step)
        listeBi.append(bi*step)
        listeBs.append(bs*step)
        listeMed1.append(med1*step)
        listeBi1.append(bi1*step)
        listeBs1.append(bs1*step)
        listeMed2.append(med2*step)
        listeBi2.append(bi2*step)
        listeBs2.append(bs2*step)


    xAbsciss = [j+1 for j in range(999)]
    k = 0
    while len(xAbsciss) < len(listeMed):
        xAbsciss.append(1000 + 2*k)
        k+=1
    plt.plot(xAbsciss,listeMed)
    plt.plot(xAbsciss,listeBi)
    plt.plot(xAbsciss,listeBs)
    plt.savefig("images/copymeasureOfPrevalenceBiDim.pdf")
    plt.show()
    plt.clf()
    plt.plot(xAbsciss,listeMed1)
    plt.plot(xAbsciss,listeBi1)
    plt.plot(xAbsciss,listeBs1)
    plt.savefig("images/copymeasureOfPrevalenceBiDim1.pdf")
    plt.show()
    plt.clf()
    plt.plot(xAbsciss,listeMed2)
    plt.plot(xAbsciss,listeBi2)
    plt.plot(xAbsciss,listeBs2)
    plt.savefig("images/copymeasureOfPrevalenceBiDim2.pdf")
    plt.show()
    plt.clf()
    
    listeMedCrop = [[],[]]
    listeMed1Crop = [[],[]]
    listeMed2Crop = [[],[]]
    for j in range(len(listeMed)):
        if listeMed[j]<0.1:
            listeMedCrop[0].append(xAbsciss[j])
            listeMedCrop[1].append(listeMed[j])
        if listeMed1[j]<0.1:
            listeMed1Crop[0].append(xAbsciss[j])
            listeMed1Crop[1].append(listeMed1[j])
        if listeMed2[j]<0.1:
            listeMed2Crop[0].append(xAbsciss[j])
            listeMed2Crop[1].append(listeMed2[j])
    plt.plot(listeMedCrop[0],listeMedCrop[1])
    plt.plot(listeMed1Crop[0],listeMed1Crop[1])
    plt.plot(listeMed2Crop[0],listeMed2Crop[1])
    plt.savefig("images/copymeasurePrevalenceCropped.pdf")
    plt.show()
    plt.clf()

    plt.plot(xAbsciss,listeMed)
    plt.plot(xAbsciss,listeBi)
    plt.plot(xAbsciss,listeBs)
    plt.plot(xAbsciss,listeMed1)
    plt.plot(xAbsciss,listeBi1)
    plt.plot(xAbsciss,listeBs1)
    plt.plot(xAbsciss,listeMed2)
    plt.plot(xAbsciss,listeBi2)
    plt.plot(xAbsciss,listeBs2)
    plt.savefig("images/copymeasureOfPrevalenceBiDimTot.pdf")
    plt.show()
    plt.clf()
    
    listeI = [listeBs[j] - listeBi[j] for j in range(len(listeBs))]
    listeI1 = [listeBs1[j] - listeBi1[j] for j in range(len(listeBs))]
    listeI2 = [listeBs2[j] - listeBi2[j] for j in range(len(listeBs))]
    plt.semilogy(xAbsciss,listeI)
    plt.semilogy(xAbsciss,listeI1)
    plt.semilogy(xAbsciss,listeI2)
    plt.savefig("images/copywidthOfCIsemilog.pdf")
    plt.show()
    plt.clf()

    fichier = open('csvForImage/measureOfPrevalenceBiDim.csv','w')
    fichier.write("Evolution of the estimation of a population containing two subpopulations :\n Population 1 has proportion {} and prevalence {}\n Population 2 has proportion {} and prevalence {}\n  Maximal group size is {}".format(ratio, prev1, 1-ratio, prev2,Nmax))
    fichier.write("Total pop :\n")
    fichier.write("Median : ")
    for el in listeMed:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Inf : ")
    for el in listeBi:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Sup : ")
    for el in listeBs:
        fichier.write(str(el)+" , ")
    fichier.write('\n')
    fichier.write("Pop1 :\n")
    fichier.write("Median : ")
    for el in listeMed1:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Inf : ")
    for el in listeBi1:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Sup : ")
    for el in listeBs1:
        fichier.write(str(el)+" , ")
    fichier.write('\n')
    fichier.write("Pop 2 :\n")
    fichier.write("Median : ")
    for el in listeMed2:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Inf : ")
    for el in listeBi2:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Sup : ")
    for el in listeBs2:
        fichier.write(str(el)+" , ")
    fichier.write('\n')
    fichier.close()

    return listeMed, listeBi, listeBs, listeMed1, listeBi1, listeBs1, listeMed2, listeBi2, listeBs2


def multinorm(listeValeurs, listeSupport, ratio, step):
    s = sum(listeValeurs)
    l = len(listeValeurs)
    d = len(listeSupport[0])
    nbPts = int(floor(1/step))+1
    
    listesProjetees = [[0 for j in range(nbPts)] for k in range(d+1)]
    
    
    for j in range(l):
        listeValeurs[j]=listeValeurs[j]/s
        ave = 0
        for i in range(d):
            k = int(floor(listeSupport[j][i]/step))
            listesProjetees[i][k] +=listeValeurs[j]
            ave += listeSupport[j][i]*ratio[i]
        k = int(floor(ave/step))
        listesProjetees[d][k]+= listeValeurs[j]

    listeMeds=[]
    listeBis=[]
    listeBss = []
    for i in range(d+1):
        med, bs, bi = normalize(listesProjetees[i])
        listeMeds.append(med*step)
        listeBss.append(bs*step)
        listeBis.append(bi*step)
    return listeMeds, listeBss, listeBis, listesProjetees

def choiceSubPop(proportions,N):
    d = len(proportions)
    result = [0 for j in range(d)]
    for k in range(N):
        s = 0
        u = npr.rand()
        v = -1
        while u > s:
            v +=1
            s += proportions[v]
        result[v] +=1
    return result
    
    
def measureOfPrevalenceMultidim(prevalences,proportions, nombreTests,Nmax,numberPoints,step):
    """    
    Simulation of a Bayesian measure of prevalence in an inhomogeneous
    population using Group testing. We optimize the group size according
    to the current knowledge of the overall prevalence

    Parameters
    ----------
    prevalences : float list
        Probability for an individual to be contaminated in a given subpopulation.
    proportions : float list
        Proportion of individuals of each subpopulation
    nombreTests : int
        Number of tests made during the experiment.
    Nmax : int
        Maximal number of individuals to sample in a single group
    numberPoints : int
        Number of sample point followed by the multidimensionalPrior
    step : float
        Width of the grid on which the prior is taken for each element

    Returns
    -------
    listeMed : float list
        Median of the prior as a function of number of tests made
    listeBi : float list
        5% quantile of the prior as a function of number of tests made.
    listeBs : float list
        95% quantile of the prior as a function of number of tests made.

    """
    numberSubPop = len(proportions)
    prevalence = 0
    for k in range(numberSubPop):
        prevalence += proportions[k]*prevalences[k]

    loiPrior = []
    supportPrior = []
    for j in range(numberPoints):
        supportPrior.append(npr.rand(numberSubPop))
        startProb=1
        for el in supportPrior[j]:
            startProb= startProb*el*(1-el)
        loiPrior.append(startProb)
    listeMed, listeBs, listeBi, listesProjetees = multinorm(loiPrior,supportPrior,proportions,step)
    
    resMed = []
    resBi = []
    resBs = []
    for d in range(numberSubPop+1):
        resMed.append([listeMed[d]])
        resBi.append([listeBi[d]])
        resBs.append([listeBs[d]])
    
    listetests = [j+1 for j in range(nombreTests)]

    for test in listetests:
        print(test)
        N = chooseSize(listesProjetees[numberSubPop],step)
        N = min(N,Nmax)

        
        representation = choiceSubPop(proportions,N)
        probNegTest = 1
        for d in range(numberSubPop):
            probNegTest = probNegTest*((1-prevalences[d])**representation[d])
        resTest = 1*(npr.rand() > probNegTest)
        
        for k in range(numberPoints):
            eps = 1
            theta = supportPrior[k]
            for pop in range(numberSubPop):
                eps = eps*(1 - theta[pop])**representation[pop]
            if resTest == 1:
                loiPrior[k] = loiPrior[k]*(1-eps)
            else:
                loiPrior[k] = loiPrior[k]*eps
        lMed, lBs, lBi, listesProjetees = multinorm(loiPrior,supportPrior,proportions,step)
        
        for pop in range(numberSubPop+1):
            resMed[pop].append(lMed[pop])
            resBi[pop].append(lBi[pop])
            resBs[pop].append(lBs[pop])

    for d in range(numberSubPop+1):
        plt.plot(resMed[d])
    plt.savefig("images/copymeasureOfPrevalenceMultiDimMedians.pdf")
    plt.show()
    plt.clf()

    for d in range(numberSubPop+1):
        plt.plot(resMed[d])
        plt.plot(resBi[d])
        plt.plot(resBs[d])
    plt.savefig("images/copymeasureOfPrevalenceMultiDim.pdf")
    plt.show()
    plt.clf()

    fichier = open('csvForImage/measureOfPrevalenceMultiDim.csv','w')
    fichier.write("Evolution of the estimation of the prevalence {} with maximal group size {}".format(prevalence,Nmax))
    fichier.write("Median : ")
    for el in resMed[numberSubPop]:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Inf : ")
    for el in resBi[numberSubPop]:
        fichier.write(str(el)+" , ")
    fichier.write("\nBorne Sup : ")
    for el in resBs[numberSubPop]:
        fichier.write(str(el)+" , ")
    fichier.write('\n\n')
    for d in range(numberSubPop):
        fichier.write("Population {} with prevalence {} and proportion {}\n".format(d+1,prevalences[d],proportions[d]))
        for el in resMed[d]:
            fichier.write(str(el)+" , ")
        fichier.write("\nBorne Inf : ")
        for el in resBi[d]:
            fichier.write(str(el)+" , ")
        fichier.write("\nBorne Sup : ")
        for el in resBs[d]:
            fichier.write(str(el)+" , ")
        fichier.write('\n\n')
    fichier.close()           

    return resMed, resBi, resBs






# listeMed, listeBi, listeBs, listeMed1, listeBi1, listeBs1, listeMed2, listeBi2, listeBs2 = measureOfPrevalenceBidim(0.005,0.05,.8,2000,200,0.0005)

measureOfPrevalenceBidim(0.005,0.05,.8,1500,200,0.0005)

# measureOfPrevalenceMultidim([0.01,0.02,0.03],[0.3,0.4,0.3],2000,200,10000,0.01)



# for prev in [0.01,0.03,0.05,0.15]:
#     print(prev)
#     listeMed, listeBi, listeBs = measureOfPrevalenceUnidim(prev,2000,50,0.00001)
#     # plt.xlabel("Number of tests")
#     # plt.ylabel("Width of credibility interval")
#     plt.semilogy([listeBs[k]-listeBi[k] for k in range(len(listeBs))])
# plt.savefig("images/copyimgCovidCommon.pdf")
# plt.show()
# plt.clf()

