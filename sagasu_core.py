#!/dls_sw/i23/scripts/ctrl_conda python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 14:48:45 2020
@author: Christian M. Orr
"""

from datetime import datetime
import os
import re
import pandas as pd
from matplotlib.patches import Ellipse
from sklearn.mixture import BayesianGaussianMixture as bgm
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import pickle
import glob

# from cuml import DBSCAN as cumlDBSCAN
from sklearn.cluster import DBSCAN
import heapq
from mpl_toolkits.mplot3d import Axes3D
import shutil
from pathlib import Path
import subprocess
import time
from tqdm import tqdm
from multiprocessing import Pool
import drmaa2
from drmaa2 import JobSession, JobTemplate, JobInfo, Drmaa2Exception


sns.set()


class core:
    def __init__(self):
        self.timestamp = datetime.now()
        self.path = os.getcwd()

    def get_input(self):
        self.projname = input("Name of project (SHELX prefix): ")
        self.fa_path = input("Path to SHELXC outputs: ")
        self.highres = int(10 * float(input("High resolution cutoff for grid: ")))
        self.lowres = int(10 * float(input("Low resolution cutoff for grid: ")))
        self.highsites = int(input("Maximum number of sites to search: "))
        self.lowsites = int(input("Minimum number of sites to search: "))
        self.ntry = int(input("Number of trials: "))
        self.clust = str(input("Run on (c)luster or (l)ocal machine? c/l ")).lower()
        self.clusteranalysis = str(
            input("Run cluster analysis after (time consuming)? y/n ")
        ).lower()
        self.insin = os.path.join(self.fa_path, self.projname + "_fa.ins")
        self.hklin = os.path.join(self.fa_path, self.projname + "_fa.hkl")
        self.writepickle()
        return (
            self.projname,
            self.fa_path,
            self.highres,
            self.lowres,
            self.highsites,
            self.lowsites,
            self.ntry,
        )

    def drmaa2template(self, workpath, jobname="sagasu"):
        # envs = {
        #     'shelx1': '_LMFILES_=/dls_sw/apps/Modules/modulefiles/global/directories:/dls_sw/apps/Modules/modulefiles/R/3.2.2:/dls_sw/apps/Modules/modulefiles/ccp4/7.1.015:/dls_sw/apps/Modules/modulefiles/shelx/ccp4',
        #     'shelx2': 'LOADEDMODULES=global/directories:R/3.2.2:ccp4/7.1.015:shelx/ccp4',
        # }
        jt = JobTemplate(
            {
                "job_name": jobname,
                "job_category": "i23_chris",
                "remote_command": "/dls/science/groups/i23/scripts/chris/Sagasu/shelxd.sh",
                "args": [str(self.projname + "_fa")],
                "min_slots": 20,
                "max_slots": 40,
                # "job_environment": envs,
                "working_directory": str(workpath),
                "output_path": str(workpath),
                "error_path": str(workpath),
                "queue_name": "low.q",
                "implementation_specific": {"uge_jt_pe": "smp",},
            }
        )
        return jt

    def drmaa2_check(self):
        job_list = [job_info[0] for job_info in self.job_details]
        time.sleep(10)
        self.session.wait_all_started(job_list)
        time.sleep(10)
        self.session.wait_all_terminated(job_list)
        
    def writepickle(self):
        with open("inps.pkl", "wb") as f:
            pickle.dump(
                [
                    self.projname,
                    self.lowres,
                    self.highres,
                    self.lowsites,
                    self.highsites,
                    self.ntry,
                    self.clusteranalysis,
                    self.clust,
                    self.insin,
                    self.hklin,
                ],
                f,
            )

    def readpickle(self):
        with open("inps.pkl", "rb") as f:
            (
                [
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
                ]
            ) = pickle.load(f)
            (
                self.projname,
                self.lowres,
                self.highres,
                self.lowsites,
                self.highsites,
                self.ntry,
                self.clusteranalysis,
                self.clust,
                self.insin,
                self.hklin,
            ) = (
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
            )
        return (
            self.projname,
            self.lowres,
            self.highres,
            self.lowsites,
            self.highsites,
            self.ntry,
            self.clusteranalysis,
            self.clust,
            self.insin,
            self.hklin,
        )

    def replace(self, file, pattern, subst):
        file_handle = open(file, "r")
        file_string = file_handle.read()
        file_handle.close()
        file_string = re.sub(pattern, subst, file_string)
        file_handle = open(file, "w")
        file_handle.write(file_string)
        file_handle.close()

    # def shelx_write(self):
    #     shelxjob = open("shelxd_job.sh", "w")
    #     shelxjob.write("module load shelx\n")
    #     shelxjob.write("shelxd " + self.projname + "_fa")
    #     shelxjob.close()
    #     os.chmod("shelxd_job.sh", 0o775)

    def run_sagasu_proc(self):
        self.session = JobSession()
        os.chdir(self.path)
        self.job_details = []
        Path(self.projname).mkdir(parents=True, exist_ok=True)
        i = self.highres
        if self.clust == "l":
            tot = (self.lowres - self.highres) * ((self.highsites + 1) - self.lowsites)
            pbar = tqdm(desc="SHELXD", total=tot, dynamic_ncols=True)
        else:
            pass
        while not (i >= self.lowres):
            Path(os.path.join(self.projname, str(i))).mkdir(parents=True, exist_ok=True)
            i2 = i / 10
            j = self.highsites
            while not (j <= (self.lowsites - 1)):
                os.makedirs(os.path.join(self.projname, str(i), str(j)), exist_ok=True)
                shutil.copy2(self.insin, (os.path.join(self.projname, str(i), str(j))))
                shutil.copy2(self.hklin, (os.path.join(self.projname, str(i), str(j))))
                workpath = os.path.join(self.path, self.projname, str(i), str(j))
                f = os.path.join(
                    self.path, self.projname, str(i), str(j), self.projname + "_fa.ins"
                )
                self.replace(f, "FIND", "FIND " + str(j) + "\n")
                self.replace(f, "SHEL", "SHEL 999 " + str(i2) + "\n")
                self.replace(f, "NTRY", "NTRY " + str(self.ntry) + "\n")
                if self.clust == "l":
                    os.chdir(workpath)
                    subprocess.run(
                        ["shelxd", self.projname + "_fa"], stdout=subprocess.PIPE
                    )
                    pbar.update(1)
                    pbar.refresh()
                    os.chdir(self.path)
                elif self.clust == "c":
                    template = self.drmaa2template(workpath)
                    job = self.session.run_job(template)
                    self.job_details.append([job])
                else:
                    print("error in input...")
                j = j - 1
            i = i + 1

    def cleanup_prev(self):
        resultspath = os.path.join(self.path, self.projname + "_results")
        self.torun = []
        if os.path.exists(resultspath):
            shutil.rmtree(resultspath)
        if not os.path.exists(self.projname + "_results"):
            os.mkdir(self.path + "/" + self.projname + "_results")
        figspath = os.path.join(self.path, self.projname + "_figures")
        if os.path.exists(figspath):
            shutil.rmtree(figspath)
        if not os.path.exists(self.projname + "_figures"):
            os.mkdir(self.path + "/" + self.projname + "_figures")
        i = self.highres
        while not (i >= self.lowres):
            j = self.highsites
            while not (j <= (self.lowsites - 1)):
                lstfile = os.path.join(
                    self.path,
                    self.projname
                    + "/"
                    + str(i)
                    + "/"
                    + str(j)
                    + "/"
                    + self.projname
                    + "_fa.lst",
                )
                self.torun.append((lstfile, i, j))
                # results(lstfile, self.path, self.projname, i, j)
                j = j - 1
            i = i + 1
        torun = self.torun
        return torun

    def results(self, filename, i, j):
        with open(filename, "r") as file:
            filedata = file.read()
            filedata = filedata.replace("/", " ")
            filedata = filedata.replace(",", " ")
            filedata = filedata.replace("CC", "")
            filedata = filedata.replace("All", "")
            filedata = filedata.replace("Weak", "")
            filedata = filedata.replace("CFOM", "")
            filedata = filedata.replace("best", "")
            filedata = filedata.replace("PATFOM", "")
            filedata = filedata.replace("CPU", "")
        with open(filename, "w") as file:
            file.write(filedata)
        with open(filename, "r") as infile, open(
            self.path
            + "/"
            + self.projname
            + "_results/"
            + str(i)
            + "_"
            + str(j)
            + ".csv",
            "w",
        ) as outfile:
            for line in infile:
                if line.startswith(" Try"):
                    line.replace("Try", "Try ")
                    outfile.write(",".join(line.split()) + "\n")
        with open(
            self.path
            + "/"
            + self.projname
            + "_results/"
            + str(i)
            + "_"
            + str(j)
            + ".csv",
            "r",
        ) as f:
            data = f.read()
            with open(
                self.path
                + "/"
                + self.projname
                + "_results/"
                + str(i)
                + "_"
                + str(j)
                + ".csv",
                "w",
            ) as w:
                w.write(data[:-1])

    def run_sagasu_analysis(self):
        clustering_distance_torun = []
        dbscan_torun = []
        hexplots_torun = []
        ccoutliers_torun = []
        if not os.path.exists(self.projname + "_figures"):
            os.mkdir(self.projname + "_figures")
        i = self.highres
        while not (i >= self.lowres):
            i2 = i / 10
            j = self.highsites
            while not (j <= (self.lowsites - 1)):
                csvfile = os.path.join(
                    self.path,
                    self.projname + "_results/" + str(i) + "_" + str(j) + ".csv",
                )
                numbers = str(i) + "_" + str(j)
                if self.clusteranalysis == "y":
                    clustering_distance_torun.append((csvfile, numbers))
                    dbscan_torun.append((csvfile, numbers, i, j))
                    hexplots_torun.append((csvfile, numbers))
                else:
                    print("No cluster analysis requested")
                ccoutliers_torun.append((csvfile, i, j))
                j = j - 1
            i = i + 1
        return clustering_distance_torun, dbscan_torun, hexplots_torun, ccoutliers_torun

    def for_ML_analysis(self):
        to_run_ML = []
        if not os.path.exists(self.projname + "_figures"):
            os.mkdir(self.projname + "_figures")
        i = self.highres
        while not (i >= self.lowres):
            i2 = i / 10
            j = self.highsites
            while not (j <= (self.lowsites - 1)):
                csvfile = os.path.join(
                    self.path,
                    self.projname + "_results/" + str(i) + "_" + str(j) + ".csv",
                )
                numbers = str(i) + "_" + str(j)
                to_run_ML.append((csvfile, numbers))
                j = j - 1
            i = i + 1
        return to_run_ML

    def plot_for_ML(self, filename, nums):
        df = pd.read_csv(
            filename,
            sep=",",
            names=[
                "linebeg",
                "TRY",
                "CPUNO",
                "CCALL",
                "CCWEAK",
                "CFOM",
                "BEST",
                "PATFOM",
            ],
        )
        plt.scatter(df["CCWEAK"], df["CCALL"], c=df["PATFOM"], cmap="Blues", marker="o")
        # plt.axis("off")
        plt.draw()
        ccallvsccweak = plt.gcf()
        ccallvsccweak.savefig(
            self.path
            + "/"
            + self.projname
            + "_figures/"
            + self.projname
            + "_"
            + nums
            + "_ML.png",
            dpi=500,
            bbox_inches=0,
        )
        ccallvsccweak.clear()
        plt.close(ccallvsccweak)

    def draw_ellipse(self, position, covariance, ax=None, **kwargs):
        """Draw an ellipse with a given position and covariance"""
        ax = ax or plt.gca()
        # Convert covariance to principal axes
        if covariance.shape == (2, 2):
            U, s, Vt = np.linalg.svd(covariance)
            angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
            width, height = 2 * np.sqrt(s)
        else:
            angle = 0
            width, height = 2 * np.sqrt(covariance)
        # Draw the Ellipse
        for nsig in range(1, 4):
            ax.add_patch(
                Ellipse(position, nsig * width, nsig * height, angle, **kwargs)
            )

    def plot_gmm(self, gmm, X, n_init, nums, label=True, ax=None):
        ax = ax or plt.gca()
        labels = gmm.fit(X).predict(X)
        if label:
            ax.scatter(X[:, 0], X[:, 1], c=labels, s=40, cmap="viridis", zorder=2)
        else:
            ax.scatter(X[:, 0], X[:, 1], s=40, zorder=2)
        w_factor = 0.2 / gmm.weights_.max()
        for pos, covar, w in zip(gmm.means_, gmm.covariances_, gmm.weights_):
            self.draw_ellipse(pos, covar, alpha=w * w_factor)
        if gmm.converged_ is True:
            print("Clustering converged after " + str(gmm.n_iter_) + " iterations")
        if gmm.converged_ is False:
            print("Clustering did not converge after " + str(n_init) + " iterations")
        print(
            str(round((gmm.weights_[0]) * 100, 1))
            + "% in first cluster, "
            + str(round((gmm.weights_[1]) * 100, 1))
            + "% in second cluster"
        )
        meanchange = np.vstack((gmm.means_, gmm.mean_prior_))
        dist_c = math.sqrt(
            ((abs((gmm.means_[0][0]) - (gmm.means_[1][0]))) ** 2)
            + ((abs((gmm.means_[0][1]) - (gmm.means_[1][1]))) ** 2)
        )
        print("Distance between clusters = " + str(dist_c))
        separation = open(self.projname + "_results/clusterseparations.csv", "a")
        separation.write(
            str(meanchange[0])
            + ","
            + str(meanchange[1])
            + ","
            + str(meanchange[2])
            + ","
            + str(dist_c)
            + ","
            + str(round((gmm.weights_[0]) * 100, 1))
            + ","
            + str(round((gmm.weights_[1]) * 100, 1))
            + "\n"
        )
        ax = plt.gcf()
        ax.savefig(
            self.path + "/" + self.projname + "_figures/" + nums + "_clsdst.png",
            dpi=300,
        )
        ax.clear()
        plt.close(ax)

    def clustering_distance(self, csvfile, nums):
        df = pd.read_csv(
            csvfile,
            sep=",",
            names=[
                "linebeg",
                "TRY",
                "CPUNO",
                "CCALL",
                "CCWEAK",
                "CFOM",
                "BEST",
                "PATFOM",
            ],
        )
        arr = df[["CCALL", "CCWEAK"]].to_numpy()
        cmean = arr.mean(axis=0)
        csd = arr.std(axis=0)
        outliermask = ((arr[:, 0]) > (cmean[0] - (2 * csd[0]))) & (
            (arr[:, 1]) > (cmean[1] - (2 * csd[1]))
        )
        arr_out = arr[outliermask]
        ni = 1000
        gmm = bgm(
            n_components=2,
            covariance_type="full",
            max_iter=ni,
            init_params="kmeans",
            tol=1e-6,
        )
        self.plot_gmm(gmm, arr, ni, nums)

    def analysis(self, filename, nums, a_res, a_sites):
        df = pd.read_csv(
            filename,
            sep=",",
            names=[
                "linebeg",
                "TRY",
                "CPUNO",
                "CCALL",
                "CCWEAK",
                "CFOM",
                "BEST",
                "PATFOM",
            ],
        )
        ccallweak = df[["CCALL", "CCWEAK"]]
        clustr = DBSCAN(eps=0.7, min_samples=1, n_jobs=-1).fit(ccallweak)
        labels = len(set(clustr.labels_))
        print("DBSCAN found " + str(labels) + " cluster(s)")
        plt.scatter(
            df["CCALL"],
            df["CCWEAK"],
            c=clustr.labels_.astype(float),
            marker="+",
            s=50,
            alpha=0.5,
        )
        plt.xlabel("CCALL")
        plt.ylabel("CCWEAK")
        plt.title("Resolution: " + str(a_res / 10) + "Å , Sites: " + str(a_sites))
        plt.draw()
        ccallvsccweak = plt.gcf()
        ccallvsccweak.savefig(
            self.path + "/" + self.projname + "_figures/" + nums + ".png", dpi=300
        )
        ccallvsccweak.clear()
        plt.close(ccallvsccweak)

    def analysis_2(self, filename, nums):
        df = pd.read_csv(
            filename,
            sep=",",
            names=[
                "linebeg",
                "TRY",
                "CPUNO",
                "CCALL",
                "CCWEAK",
                "CFOM",
                "BEST",
                "PATFOM",
            ],
        )
        sns.jointplot(x=df["CCALL"], y=df["CCWEAK"], kind="hex", space=0)
        plt.draw()
        snsplot = plt.gcf()
        snsplot.savefig(
            self.path + "/" + self.projname + "_figures/" + nums + "_hexplot.png",
            dpi=300,
        )
        snsplot.clear()
        plt.close(snsplot)

    def CFOM_PATFOM_analysis(self, filename, resolution, sitessearched):
        df = pd.read_csv(
            filename,
            sep=",",
            names=[
                "linebeg",
                "TRY",
                "CPUNO",
                "CCALL",
                "CCWEAK",
                "CFOM",
                "BEST",
                "PATFOM",
            ],
        )
        pd.DataFrame.drop(df, labels="linebeg", axis=1, inplace=True)
        df.sort_values("CFOM", ascending=False, inplace=True, na_position="last")
        top_CFOM = df["CFOM"].values[0]
        corr_PATFOM = df["PATFOM"].values[0]
        with open(self.projname + "_results/CFOM_PATFOM.csv", "a") as allfom:
            allfom.write(
                str(int(resolution) / 10)
                + ","
                + str(sitessearched)
                + ","
                + str(top_CFOM)
                + ","
                + str(corr_PATFOM)
                + "\n"
            )

    def ccalloutliers(self, filename, resolution, sitessearched):
        df = pd.read_csv(
            filename,
            sep=",",
            names=[
                "linebeg",
                "TRY",
                "CPUNO",
                "CCALL",
                "CCWEAK",
                "CFOM",
                "BEST",
                "PATFOM",
            ],
        )
        pd.DataFrame.drop(df, labels="linebeg", axis=1, inplace=True)
        median = df["CCALL"].median()
        arr = df[["CCALL", "CCWEAK"]].to_numpy()
        cmean = arr.mean(axis=1)
        csd = arr.std(axis=1)
        outliermask = ((arr[:, 0]) > (cmean[0] - (2 * csd[0]))) & (
            (arr[:, 1]) > (cmean[1] - (2 * csd[1]))
        )
        arr = arr[outliermask]
        mad = np.median(np.sqrt((arr[:, 0] - median) ** 2))
        ccallmad = arr[:, 0] - median
        mad10 = sum(i > 10 * mad for i in ccallmad)
        mad9 = sum(i > 9 * mad for i in ccallmad)
        mad8 = sum(i > 8 * mad for i in ccallmad)
        mad7 = sum(i > 7 * mad for i in ccallmad)
        mad6 = sum(i > 6 * mad for i in ccallmad)
        mad5 = sum(i > 5 * mad for i in ccallmad)
        allmad = open(self.projname + "_results/ccall.csv", "a")
        allmad.write(
            str(int(resolution) / 10)
            + ","
            + str(sitessearched)
            + ","
            + str(mad5)
            + ","
            + str(mad6)
            + ","
            + str(mad7)
            + ","
            + str(mad8)
            + ","
            + str(mad9)
            + ","
            + str(mad10)
            + "\n"
        )
        allmad.close()

    def ccweakoutliers(self, filename, resolution, sitessearched):
        df = pd.read_csv(
            filename,
            sep=",",
            names=[
                "linebeg",
                "TRY",
                "CPUNO",
                "CCALL",
                "CCWEAK",
                "CFOM",
                "BEST",
                "PATFOM",
            ],
        )
        pd.DataFrame.drop(df, labels="linebeg", axis=1, inplace=True)
        median = df["CCWEAK"].median()
        arr = df[["CCALL", "CCWEAK"]].to_numpy()
        cmean = arr.mean(axis=0)
        csd = arr.std(axis=0)
        outliermask = ((arr[:, 0]) > (cmean[0] - (2 * csd[0]))) & (
            (arr[:, 1]) > (cmean[1] - (2 * csd[1]))
        )
        arr = arr[outliermask]
        mad = np.median(np.sqrt((arr[:, 1] - median) ** 2))
        ccweakmad = arr[:, 1] - median
        mad10 = sum(i > 10 * mad for i in ccweakmad)
        mad9 = sum(i > 9 * mad for i in ccweakmad)
        mad8 = sum(i > 8 * mad for i in ccweakmad)
        mad7 = sum(i > 7 * mad for i in ccweakmad)
        mad6 = sum(i > 6 * mad for i in ccweakmad)
        mad5 = sum(i > 5 * mad for i in ccweakmad)
        allmad = open(self.projname + "_results/ccweak.csv", "a")
        allmad.write(
            str(int(resolution) / 10)
            + ","
            + str(sitessearched)
            + ","
            + str(mad5)
            + ","
            + str(mad6)
            + ","
            + str(mad7)
            + ","
            + str(mad8)
            + ","
            + str(mad9)
            + ","
            + str(mad10)
            + "\n"
        )
        allmad.close()

    def tophits(self):
        df = pd.read_csv(
            self.projname + "_results/ccall.csv",
            sep=",",
            names=["res", "sites", "mad5", "mad6", "mad7", "mad8", "mad9", "mad10"],
        )
        df["score"] = (
            (df["mad5"] * 1)
            + (df["mad6"] * 4)
            + (df["mad7"] * 8)
            + (df["mad8"] * 32)
            + (df["mad9"] * 128)
            + (df["mad10"] * 512)
        )
        df.sort_values(by=["score"], ascending=False, inplace=True)
        top = df[["res", "sites", "score"]]
        top = df.head(10)
        self.topallhtml = top.reset_index(drop=True).to_html()
        self.topall = str(top.reset_index(drop=True))
        weak_df = pd.read_csv(
            self.projname + "_results/ccweak.csv",
            sep=",",
            names=["res", "sites", "mad5", "mad6", "mad7", "mad8", "mad9", "mad10"],
        )
        weak_df["score"] = (
            (weak_df["mad5"] * 1)
            + (weak_df["mad6"] * 4)
            + (weak_df["mad7"] * 8)
            + (weak_df["mad8"] * 32)
            + (weak_df["mad9"] * 128)
            + (weak_df["mad10"] * 512)
        )
        weak_df.sort_values(by=["score"], ascending=False, inplace=True)
        top = weak_df[["res", "sites", "score"]]
        top = weak_df.head(10)
        self.topweakhtml = top.reset_index(drop=True).to_html()
        self.topweak = str(top.reset_index(drop=True))
        cfom_df = pd.read_csv(
            self.projname + "_results/CFOM_PATFOM.csv",
            sep=",",
            names=["res", "sites", "CFOM", "PATFOM"],
        )
        cfom_df.sort_values("CFOM", ascending=False, inplace=True)
        cfom_df["score"] = (
            (cfom_df["CFOM"])
            - (((cfom_df["res"])) * (cfom_df["res"]))
            - (0.3 * (cfom_df["sites"]))
        )
        top = cfom_df.head(10)
        self.top_CFOMhtml = top.reset_index(drop=True).to_html()
        self.top_CFOM = str(top.reset_index(drop=True))
        with open("tophits.txt", "w") as outfile:
            outfile.write(self.topall)
            outfile.write("\n")
            outfile.write(self.topweak)
            outfile.write("\n")
            outfile.write(self.top_CFOM)
        # make some 3d figures
        ax = plt.axes(projection="3d")
        ax.plot_trisurf(
            df["res"], df["sites"], df["score"], cmap="viridis", edgecolor="none"
        )
        madplot = plt.gcf()
        madplot.savefig(self.projname + "_figures/ccall.png", dpi=600)
        plt.clf()
        plt.cla()
        plt.close()
        ax = plt.axes(projection="3d")
        ax.plot_trisurf(
            weak_df["res"],
            weak_df["sites"],
            weak_df["score"],
            cmap="viridis",
            edgecolor="none",
        )
        madplot = plt.gcf()
        madplot.savefig(self.projname + "_figures/ccweak.png", dpi=600)
        plt.clf()
        plt.cla()
        plt.close()
        ax = plt.axes(projection="3d")
        ax.plot_trisurf(
            cfom_df["res"],
            cfom_df["sites"],
            cfom_df["score"],
            cmap="viridis",
            edgecolor="none",
        )
        madplot = plt.gcf()
        madplot.savefig(self.projname + "_figures/CFOM.png", dpi=600)
        # run phenix.emma on top 2 CCALL
        top = df[["res", "sites", "score"]]
        top = df.head(10)
        (firstres, firstsites, secondres, secondsites) = (
            top.iloc[[0], [0]].values[0],
            top.iloc[[0], [1]].values[0],
            top.iloc[[1], [0]].values[0],
            top.iloc[[1], [1]].values[0],
        )
        (firstres, firstsites, secondres, secondsites) = (
            ((firstres * 10).astype(np.int)).item(0),
            (firstsites.astype(np.int)).item(0),
            ((secondres * 10).astype(np.int)).item(0),
            (secondsites.astype(np.int)).item(0),
        )
        with open(
            self.path
            + "/"
            + self.projname
            + "/"
            + str(firstres)
            + "/"
            + str(firstsites)
            + "/"
            + self.projname
            + "_fa.res",
            "r",
        ) as infile:
            for line in infile:
                if line.startswith("TITL"):
                    words = line.split()
                    self.sg = words[-1]
        self.emma = os.popen(
            "phenix.emma "
            + str(
                os.path.join(
                    self.path,
                    self.projname,
                    str(firstres),
                    str(firstsites),
                    (self.projname + "_fa.pdb "),
                )
                + os.path.join(
                    self.path,
                    self.projname,
                    str(secondres),
                    str(secondsites),
                    (self.projname + "_fa.pdb"),
                )
                + " --tolerance=6 --space_group="
                + self.sg
            )
        ).read()
        self.emmain = str(
            "First model - "
            + str(float(firstres / 10))
            + " Å with a sites cutoff of "
            + str(firstsites)
            + "\n"
            + "Second model - "
            + str(float(secondres / 10))
            + " Å with a sites cutoff of "
            + str(secondsites)
        )

    def writehtml(self):
        self.html_init = """
        <!doctype html>
        <html>
        <head>
            <title>Sagasu - SHELXD Grid</title>
            <style>
            background: linear-gradient(0deg, #e4fffd 0%, #60f3ff 100%);
            </style>
        </head>
        <body>
        <h1 style="text-align: center;"><span style="font-family:courier new,courier,monospace;">Sagasu - SHELXD Grid Search </span></h1>

        <p><span style="font-family:courier new,courier,monospace;">Results for project <strong>{projname}</strong>, <strong>{ntry}</strong> trys with a low resolution limit of <strong>{lowres}</strong> and a high resolution limit of <strong>{highres}</strong>, searching for a number of sites between <strong>{lowsites}</strong> and <strong>{highsites}</strong>.</span></p>

        <hr />
        """.format(
            projname=self.projname,
            ntry=str(self.ntry),
            lowres=str(float(self.lowres / 10)),
            highres=str(float(self.highres / 10)),
            lowsites=str(self.lowsites),
            highsites=str(self.highsites),
        )

        self.html_topten = """
        <p><span style="font-family:courier new,courier,monospace;"><span style="font-size:18px;"><strong><u>Here are the top 10 hits:</u></strong></span></span></p>

        <p><span style="font-family:courier new,courier,monospace;"><strong>For CCALL:</strong></span></p>

        <p><span style="font-family:courier new,courier,monospace;">{CCALL_tophits}</span></p>

        <p><span style="font-family:courier new,courier,monospace;"><strong>For CCWEAK:</strong></span></p>

        <p><span style="font-family:courier new,courier,monospace;">{CCWEAK_tophits}</span></p>

        <p><span style="font-family:courier new,courier,monospace;"><strong>For CFOM:</strong></span></p>

        <p><span style="font-family:courier new,courier,monospace;">{CFOM_tophits}</span></p>

        <hr />
        <p><span style="font-family:courier new,courier,monospace;">phenix.emma output:</span></p>
        <p><span style="font-family:courier new,courier,monospace; white-space: pre-line">
        """.format(
            CCALL_tophits=str(self.topallhtml),
            CCWEAK_tophits=str(self.topweakhtml),
            CFOM_tophits=str(self.top_CFOMhtml),
        )

        with open("sagasu.html", "w") as htmlfile:
            htmlfile.write(self.html_init + "\n")
        with open("sagasu.html", "a") as htmlfile:
            htmlfile.write(self.html_topten + "\n")
            htmlfile.write(self.emmain + "\n")
            for line in self.emma.splitlines():
                htmlfile.write(line + "\n")
            htmlfile.write("</span></p>" + "\n")
            htmlfile.write(
                """<hr />
                <p><span style="font-family:courier new,courier,monospace;"><span style="font-size:18px;"><strong><u>Plots:</u></strong></span></span></p>"""
            )
            htmlfile.write(
                """
                <table><tbody><tr>
                <th><p><img title="{projname} CCALL" src="{projname}_figures/ccall.png" style="float: left; border-width: 2px; border-style: solid; width: 768px; height: 576px;" /></th>
            """.format(
                    projname=self.projname
                )
            )
            htmlfile.write(
                """
                <th><p><img title="{projname} CCWEAK" src="{projname}_figures/ccweak.png" style="float: left; border-width: 2px; border-style: solid; width: 768px; height: 576px;" /></th>
            """.format(
                    projname=self.projname
                )
            )
            htmlfile.write(
                """
                <th><p><img title="{projname} CFOM" src="{projname}_figures/CFOM.png" style="float: left; border-width: 2px; border-style: solid; width: 768px; height: 576px;" /></th>
                </tr></tbody></table>
            """.format(
                    projname=self.projname
                )
            )
            for plot in glob.glob(
                os.path.join(self.path, (self.projname + "_figures"), "*ML*")
            ):
                htmlfile.write(
                    """
                <p><img title="{plot}" src="{plot}" style="float: left; border-width: 2px; border-style: solid; width: 420px; height: 320px;" />
                """.format(
                        plot=str(plot)
                    )
                )
            htmlfile.write(
                """</p>
        </body>
        </html>"""
            )


if __name__ == "__main__":
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
        run = core()
        run.get_input()
        run.writepickle()
        if os.path.exists(os.path.join(path, "inps.pkl")):
            run.readpickle()
            # run.shelx_write(projname)
            run.run_sagasu_proc()
        if clust == "c":
            run.qstat_progress(lowres, highres, lowsites, highsites)
        else:
            print("Processing finished.")

    if pro_or_ana == "a" or "p":
        run = core()
        print("Analysis running, prepping files...")
        if os.path.exists(os.path.join(path, "inps.pkl")):
            run.readpickle()
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
            print("Generating hexplots...")
            pool.starmap(run.analysis_2, hexplots_torun)
            print("Running outlier analysis...")
            pool.starmap(run.ccalloutliers, ccoutliers_torun)
            pool.starmap(run.ccweakoutliers, ccoutliers_torun)
            run.tophits()
            to_run_ML = run.for_ML_analysis()
            pool.starmap(run.plot_for_ML, to_run_ML)
        else:
            print("No previous run found")
