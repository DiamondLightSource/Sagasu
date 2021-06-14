import sagasu_core
import os
import sys
from multiprocessing import Pool

pool = Pool(os.cpu_count() - 1)
print("Using ", str(os.cpu_count() - 1), "CPU cores")

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
) = sagasu_core.core.readpickle()

ml_plots = sagasu_core.core()
ml_plots.readpickle()

if os.path.exists(os.path.join(path, projname)):
    print("Data appears to be here, continuing...")
else:
    sys.exit()

if os.path.exists(os.path.join(path, projname + "_results")):
    to_run_ML = ml_plots.for_ML_analysis()
    pool.map(ml_plots.plot_for_ML, to_run_ML)
else:
    to_run = ml_plots.cleanup_prev()
    pool.map(ml_plots.results, to_run)
    to_run_ML = ml_plots.for_ML_analysis()
    pool.map(ml_plots.plot_for_ML, to_run_ML)
