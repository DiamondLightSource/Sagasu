#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Chris
"""
import os
import sagasu_core
import pickle
import gridmap as gm

pro_or_ana = str(input("Would you like to run (p)rocessing and analysis or just (a)nalysis: ").lower())
path = os.getcwd()

if pro_or_ana == 'p':
    projname = input("Name of project (SHELX prefix): ")
    fa_path = input("Path to SHELXC outputs: ")
    highres = input("High resolution cutoff for grid: ")
    lowres = input("Low resolution cutoff for grid: ")
    highsites = input("Maximum number of sites to search: ")
    lowsites = input("Minimum number of sites to search: ")
    ntry = input("Number of trials: ")
    clust = str(input("Run on (c)luster or (l)ocal machine? ")).lower()
    clusteranalysis = "y"
    pro_or_ana = str(pro_or_ana).lower()
    highres = int((10 * float(highres)))
    lowres = int((10 * float(lowres)))
    highsites = int(highsites)
    lowsites = int(lowsites)
    insin = os.path.join(fa_path, projname + "_fa.ins")
    hklin = os.path.join(fa_path, projname + "_fa.hkl")
    ntry = int(ntry)
    statusofrun = "-hold_jid "
    clust = str(clust).lower()

os.chdir(path)

if pro_or_ana == 'p':
    print("Processing started")
    sagasu_core.writepickle(path,
                            projname,
                            lowres,
                            highres,
                            lowsites,
                            highsites,
                            ntry,
                            clusteranalysis
                            )
    sagasu_core.shelx_write(projname)
    sagasu_core.run_sagasu_proc(pro_or_ana,
                                projname,
                                highres,
                                lowres,
                                highsites,
                                lowsites,
                                insin,
                                hklin,
                                path,
                                ntry,
                                statusofrun,
                                clust
                                )

if pro_or_ana == 'a' | 'p':
    print("Analysis running...")
    with open("inps.pkl", "rb") as f:
                (path,
                projname,
                lowres,
                highres,
                lowsites,
                highsites,
                ntry,
                clusteranalysis) = pickle.load(f)
    sagasu_core.cleanup_prev(
        path, projname, highres, lowres, highsites, lowsites
        )
    sagasu_core.run_sagasu_analysis(
        projname, highres, lowres, highsites, lowsites, path, clusteranalysis
        )
    sagasu_core.tophits(projname, path)
