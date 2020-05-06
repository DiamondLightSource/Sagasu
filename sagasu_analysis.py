#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:24:17 2020

@author: vwg85559
"""

import pickle
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
import heapq
from mpl_toolkits.mplot3d import Axes3D

with open('inps.pkl', 'rb') as f:
    path, projname, lowres, highres, lowsites, highsites, ntry, clusteranalysis = pickle.load(f)
 
def results(filename):
    with open(filename, 'r') as file:
        filedata = file.read()
        filedata = filedata.replace('/', ' ')
        filedata = filedata.replace(',', ' ')
        filedata = filedata.replace('CC', '')
        filedata = filedata.replace('All', '')
        filedata = filedata.replace('Weak', '')
        filedata = filedata.replace('CFOM', '')
        filedata = filedata.replace('best', '')
        filedata = filedata.replace('PATFOM', '')
        filedata = filedata.replace('CPU', '')
    with open(filename, 'w') as file:
        file.write(filedata)
    with open(filename, 'r') as infile, open(path+"/"+projname+"_results/"+str(i)+"_"+str(j)+".csv", 'w') as outfile:
        for line in infile:
            if line.startswith(" Try"):
                outfile.write(','.join(line.split())+'\n')
    with open(path+"/"+projname+"_results/"+str(i)+"_"+str(j)+".csv", 'r') as f:
        data = f.read()
        with open(path+"/"+projname+"_results/"+str(i)+"_"+str(j)+".csv", 'w') as w:
            w.write(data[:-1])

def analysis(filename, nums, a_res, a_sites):
    df = pd.read_csv(filename, sep=',', names=['linebeg', 'TRY', 'CPUNO', 'CCALL', 'CCWEAK', 'CFOM', 'BEST', 'PATFOM'])
    ccallweak = df[['CCALL', 'CCWEAK']]
    clustr = DBSCAN(eps=0.7, min_samples=1, n_jobs=-1).fit(ccallweak)
    labels = len(set(clustr.labels_))
    print("There are "+str(labels)+" cluster(s)")
    plt.scatter(df['CCALL'], df['CCWEAK'], c= clustr.labels_.astype(float), marker='+', s=50, alpha=0.5)
    plt.xlabel("CCALL")
    plt.ylabel("CCWEAK")
    plt.title("Resolution: "+str(a_res/10)+"Å , Sites: "+str(a_sites))
    plt.draw()
    ccallvsccweak = plt.gcf()
    ccallvsccweak.savefig(path+"/"+projname+"_figures/"+nums+".png", dpi=300)
    #plt.show()

def outliers(filename, resolution, sitessearched):
    df = pd.read_csv(filename, sep=',', names=['linebeg', 'TRY', 'CPUNO', 'CCALL', 'CCWEAK', 'CFOM', 'BEST', 'PATFOM'])
    pd.DataFrame.drop(df, labels = "linebeg", axis=1, inplace=True)
    median = df['CCALL'].median()
    mad = np.median(np.sqrt((df['CCALL'] - median)**2))
    ccallmax = heapq.nlargest(3, df['CCALL'])
    ccallmad = df['CCALL'] - median
    mad10 = sum(i > 10 * mad for i in ccallmad)
    mad9 = sum(i > 9 * mad for i in ccallmad)
    mad8 = sum(i > 8 * mad for i in ccallmad)
    mad7 = sum(i > 7 * mad for i in ccallmad)
    mad6 = sum(i > 6 * mad for i in ccallmad)
    mad5 = sum(i > 5 * mad for i in ccallmad)
    print("For a resolution cutoff of "+str(resolution)+" with "+str(sitessearched)+" sites:")
    print("number of CCall with CCall - median > 10 * MAD = "+str(mad10))
    print("number of CCall with CCall - median >  9 * MAD = "+str(mad9))
    print("number of CCall with CCall - median >  8 * MAD = "+str(mad8))
    print("number of CCall with CCall - median >  7 * MAD = "+str(mad7))
    print("number of CCall with CCall - median >  6 * MAD = "+str(mad6))
    print("number of CCall with CCall - median >  5 * MAD = "+str(mad5))
    print("Three largest CCall values: "+str(ccallmax))
    print("")
    allmad = open(projname+"_results/mad.csv", 'a') 
    allmad.write(str(int(resolution) / 10)+","+str(sitessearched)+","+str(mad5)+","+str(mad6)+","+str(mad7)+","+str(mad8)+","+str(mad9)+","+str(mad10)+"\n")
    allmad.close()   
    
################ ANALYSIS ################

#prep files from shelxd runs
if not os.path.exists(projname+"_results"):
    os.system("mkdir "+projname+"_results")
i = highres
while not (i >= lowres):
    i2 = (i/10)
    j = highsites
    while not (j <= (lowsites - 1)):
        lstfile = os.path.join(path, projname+"/"+str(i)+"/"+str(j)+"/"+projname+"_fa.lst")
        results(lstfile)
        j = j - 1
    i = i + 1 

#run cluster outlier analysis with optional cluster analysis
if not os.path.exists(projname+"_figures"):
    os.system("mkdir "+projname+"_figures")
i = highres
while not (i >= lowres):
    i2 = (i/10)
    j = highsites
    while not (j <= (lowsites - 1)):
        csvfile = os.path.join(path, projname+"_results/"+str(i)+"_"+str(j)+".csv")
        numbers = str(i)+"_"+str(j)
        if clusteranalysis == 'y':
            analysis(csvfile, numbers, i, j)
        else:
            print("No cluster analysis requested")
        outliers(csvfile, i, j)
        j = j - 1
    i = i + 1
   
#create and save top hits and madplot
df = pd.read_csv(projname+"_results/mad.csv", sep=',', names=['res', 'sites', 'mad5', 'mad6', 'mad7', 'mad8', 'mad9', 'mad10'])
df['score'] = (df['mad5'] * 1) + (df['mad6'] * 4) + (df['mad7'] * 8) + (df['mad8'] * 32) + (df['mad9'] * 128) + (df['mad10'] * 512)
df.sort_values(by=['score'], ascending=False,inplace=True)
top = df[['res', 'sites', 'score']]
top = df.head(10)
print("""
Here are the top 10 hits:

""")
print(top)
ax = plt.axes(projection='3d')
ax.plot_trisurf(df['res'], df['sites'], df['score'], cmap='viridis', edgecolor='none')
madplot = plt.gcf()
madplot.savefig(projname+"_figures/mad.png", dpi=600)