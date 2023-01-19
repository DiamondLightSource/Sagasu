#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Chris
"""
import sagasu_core
import os
from multiprocessing import Pool
from halo import Halo


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
        with Halo(
            text="\nSubmitting jobs",
            text_color="green",
            spinner="monkey",
        ):
            run.run_sagasu_proc()
        with Halo(
            text="\nJobs are running, please be patient and watch the shark",
            text_color="green",
            spinner="shark",
        ):
            run.drmaa2_check()
    else:
        pass

if pro_or_ana == "a" or "p":
    run = sagasu_core.core()
    if os.path.exists(os.path.join(path, "inps.pkl")):
        run.readpickle()
        to_run = run.cleanup_prev()
        with Halo(
            text="\nPulling out the important stuff", text_color="green", spinner="dots12"
        ):
            pool.starmap(run.results, to_run)
        (
            clustering_distance_torun,
            dbscan_torun,
            hexplots_torun,
            ccoutliers_torun,
        ) = run.run_sagasu_analysis()
        # print("Clustering distance analysis...")
        # pool.starmap(run.clustering_distance, clustering_distance_torun)
        # print("DBScan")
        # pool.starmap(run.analysis, dbscan_torun)
        # print("Generating hexplots...")
        # pool.starmap(run.analysis_2, hexplots_torun)
        with Halo(text="\nLooking for outliers", text_color="green", spinner="toggle"):
            pool.starmap(run.ccalloutliers, ccoutliers_torun)
            pool.starmap(run.ccweakoutliers, ccoutliers_torun)
            pool.starmap(run.CFOM_PATFOM_analysis, ccoutliers_torun)
            run.vectoroutliers()           
            run.tophits()
        with Halo(
            text="\nGenerating pretty pictures", text_color="green", spinner="pong"
        ):
            to_run_ML = run.for_ML_analysis()
            pool.starmap(run.plot_for_ML, to_run_ML)
        run.writehtml()
        print("\nRun 'firefox sagasu.html' to view results")
    else:
        print("No previous run found!")
