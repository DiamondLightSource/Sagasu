#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Chris
"""
import os
import sagasu_core
import pickle

pro_or_ana = str(
    input(
        "Would you like to run (p)rocessing and analysis or just (a)nalysis: "
    ).lower()
)
path = os.getcwd()

os.chdir(path)

if pro_or_ana == "p":
    sagasu_core.get_input()
    sagasu_core.writepickle(
        path, projname, lowres, highres, lowsites, highsites, ntry, clusteranalysis
    )
    sagasu_core.shelx_write(projname)
    sagasu_core.run_sagasu_proc(
        pro_or_ana,
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
        clust,
    )
    sagasu_core.qstat_progress(lowres, highres, lowsites, highsites)

if pro_or_ana == "a" or "p":
    print("Analysis running...")
    if os.path.exists(os.path.join(path, "inps.pkl")):
        with open("inps.pkl", "rb") as f:
            (
                path,
                projname,
                lowres,
                highres,
                lowsites,
                highsites,
                ntry,
                clusteranalysis,
            ) = pickle.load(f)
        sagasu_core.cleanup_prev(path, projname, highres, lowres, highsites, lowsites)
        sagasu_core.run_sagasu_analysis(
            projname, highres, lowres, highsites, lowsites, path, clusteranalysis
        )
        sagasu_core.tophits(projname, path)
    else:
        print("No previous run found here, are you sure you are in the correct path?")