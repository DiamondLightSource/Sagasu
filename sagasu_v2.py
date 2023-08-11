import os
import sagasu_core
from multiprocessing import Pool
from halo import Halo

pool = Pool(os.cpu_count() - 1)


def process_data(run):
    to_run, to_run_prasa = run.cleanup_prev()

    with Halo(
        text="Pulling out the important stuff", text_color="green", spinner="dots12"
    ):
        pool.starmap(run.results, to_run)
        pool.starmap(run.prasa_results, to_run_prasa)

    ccoutliers_torun = run.run_sagasu_analysis()

    with Halo(text="Looking for outliers", text_color="green", spinner="toggle"):
        pool.starmap(run.ccalloutliers, ccoutliers_torun)
        pool.starmap(run.ccweakoutliers, ccoutliers_torun)
        pool.starmap(run.CFOM_PATFOM_analysis, ccoutliers_torun)
        run.vectoroutliers()
        run.tophits()

    with Halo(text="Generating pretty pictures", text_color="green", spinner="pong"):
        to_run_ML = run.for_ML_analysis()
        pool.starmap(run.plot_for_ML, to_run_ML)

    run.writehtml()
    print("\nRun 'firefox sagasu.html' to view results")


def main():
    path = os.getcwd()
    print("You are here:", path)
    print("Using", os.cpu_count() - 1, "CPU cores")

    pro_or_ana = input(
        "Would you like to run (p)rocessing and analysis or just (a)nalysis: "
    ).lower()

    if pro_or_ana == "p":
        run = sagasu_core.core()
        run.writepickle()

        if os.path.exists(os.path.join(path, "inps.pkl")):
            run.readpickle()
            run.prasa_prep()
            run.shelxd_prep()

            with Halo(text="Prepping jobs", text_color="green", spinner="pipe"):
                run.run_sagasu_proc()

            with Halo(text="Submitting jobs", text_color="green", spinner="monkey"):
                run.submit_shelxd_job_slurm()
                run.submit_afroprasa_job_slurm()

            with Halo(
                text="Jobs are running, please be patient and watch the shark",
                text_color="green",
                spinner="shark",
            ):
                waiting = run.wait_for_slurm_jobs()
                if not waiting:
                    exit()

    if pro_or_ana == "a" or pro_or_ana == "p":
        run = sagasu_core.core()

        if os.path.exists(os.path.join(path, "inps.pkl")):
            run.readpickle()
            process_data(run)
        else:
            print("No previous run found!")


if __name__ == "__main__":
    main()
