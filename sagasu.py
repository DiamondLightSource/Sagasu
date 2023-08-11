#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
        run.prasa_prep()
        run.shelxd_prep()
        with Halo(
            text="\nPrepping jobs",
            text_color="green",
            spinner="pipe",
        ):
            run.run_sagasu_proc()
        with Halo(
            text="\nSubmitting jobs",
            text_color="green",
            spinner="monkey",
        ):
            run.get_slurm_token()
            run.submit_shelxd_job_slurm()
            run.submit_afroprasa_job_slurm()
        with Halo(
            text="\nJobs are running, please be patient and watch the shark",
            text_color="green",
            spinner="shark",
        ):
            waiting = run.wait_for_slurm_jobs()
            if not waiting:
                exit()
    else:
        pass

if pro_or_ana == "a" or "p":
    run = sagasu_core.core()
    if os.path.exists(os.path.join(path, "inps.pkl")):
        run.readpickle()
        to_run, to_run_prasa = run.cleanup_prev()
        with Halo(
            text="\nPulling out the important stuff",
            text_color="green",
            spinner="dots12",
        ):
            pool.starmap(run.results, to_run)
            pool.starmap(run.prasa_results, to_run_prasa)
        ccoutliers_torun = run.run_sagasu_analysis()
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
