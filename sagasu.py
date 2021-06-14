#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Chris
"""
import sagasu_core
import pickle
import os

path = os.getcwd()
print("You are here:", path)
pro_or_ana = str(
    input(
        "Would you like to run (p)rocessing and analysis or just (a)nalysis: "
    ).lower()
)

run = sagasu_core.core()

if pro_or_ana == "p":
    run.get_input()
    if os.path.exists(os.path.join(path, "inps.pkl")):
        with open("inps.pkl", "rb") as f:
            (
                projname,
                lowres,
                highres,
                lowsites,
                highsites,
                ntry,
                clusteranalysis,
                clust,
                insin,
                hklin,
            ) = pickle.load(f)
    run.shelx_write(projname)
    run.run_sagasu_proc(
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
        clust,
    )
    if clust == "c":
        run.qstat_progress(lowres, highres, lowsites, highsites)
    else:
        print("Processing finished.")

if pro_or_ana == "a" or "p":
    print("Analysis running...")
    if os.path.exists(os.path.join(path, "inps.pkl")):
        with open("inps.pkl", "rb") as f:
            (
                projname,
                lowres,
                highres,
                lowsites,
                highsites,
                ntry,
                clusteranalysis,
                clust,
                insin,
                hklin,
            ) = pickle.load(f)
        run.cleanup_prev(path, projname, highres, lowres, highsites, lowsites)
        run.run_sagasu_analysis(
            projname, highres, lowres, highsites, lowsites, path, clusteranalysis
        )
        run.tophits(projname, path)
    else:
        print("No previous run found here, are you sure you are in the correct path?")
