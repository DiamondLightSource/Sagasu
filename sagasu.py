#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Chris
"""
import sagasu_core
import pickle
import os

path = os.getcwd()

pro_or_ana = str(
    input(
        "Would you like to run (p)rocessing and analysis or just (a)nalysis: "
    ).lower()
)

if pro_or_ana == "p":
    sagasu_core.get_input()
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
        clust,
    )
    if clust == "c":
        sagasu_core.qstat_progress(lowres, highres, lowsites, highsites)
    else:
        print("Processing finished.")

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
