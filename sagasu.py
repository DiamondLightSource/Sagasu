#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Chris
"""
import sagasu_core
import os
from multiprocessing import Pool


path = os.getcwd()
print("You are here:", path)
pool = Pool(os.cpu_count() - 1)
print("Using ", str(os.cpu_count() - 1), "CPU cores")
pro_or_ana = str(
    input(
        "Would you like to run (p)rocessing and analysis or just (a)nalysis: "
    ).lower()
)


if pro_or_ana == "p":
    run = sagasu_core.core()
    projname, fa_path, highres, lowres, highsites, lowsites, ntry = run.get_input()
    run.writepickle()
    if os.path.exists(os.path.join(path, "inps.pkl")):
        run.readpickle()
        run.shelx_write()
        run.run_sagasu_proc()
        gonext = input("Continue? ")
    else:
        print("Processing finished.")

if pro_or_ana == "a" or "p":
    run = sagasu_core.core()
    print("Analysis running, prepping files...")
    if os.path.exists(os.path.join(path, "inps.pkl")):
        run.readpickle()
        to_run = run.cleanup_prev()
        pool.starmap(run.results, to_run)
        clustering_distance_torun, dbscan_torun, hexplots_torun, ccoutliers_torun = (
            run.run_sagasu_analysis()
        )
        print("Clustering distance analysis...")
        # pool.starmap(run.clustering_distance, clustering_distance_torun)
        print("DBScan")
        # pool.starmap(run.analysis, dbscan_torun)
        print("Generating hexplots...")
        # pool.starmap(run.analysis_2, hexplots_torun)
        print("Running outlier analysis...")
        pool.starmap(run.ccalloutliers, ccoutliers_torun)
        pool.starmap(run.ccweakoutliers, ccoutliers_torun)
        pool.starmap(run.CFOM_PATFOM_analysis, ccoutliers_torun)
        run.tophits()
        to_run_ML = run.for_ML_analysis()
        print("Making ML plots...")
        pool.starmap(run.plot_for_ML, to_run_ML)
        run.writehtml()
    else:
        print("No previous run found!")
