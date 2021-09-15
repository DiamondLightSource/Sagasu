#!/usr/bin/env python3
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
#from cuml import DBSCAN as cumlDBSCAN
from sklearn.cluster import DBSCAN
import heapq
from mpl_toolkits.mplot3d import Axes3D
import shutil
from pathlib import Path
import subprocess
import time
from tqdm import tqdm
from multiprocessing import Pool


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
        self.clusteranalysis = str(input("Run cluster analysis after (time consuming)? y/n ")).lower()
        self.insin = os.path.join(self.fa_path, self.projname + "_fa.ins")
        self.hklin = os.path.join(self.fa_path, self.projname + "_fa.hkl")
        self.writepickle()

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
            self.projname, self.lowres, self.highres, self.lowsites, self.highsites, self.ntry, self.clusteranalysis, self.clust, self.insin, self.hklin = (
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

    def qstat_progress(self):
        runs = (((self.lowres - 1) - self.highres)) * (self.highsites - self.lowsites)
        q = subprocess.Popen("qstat", stdout=subprocess.PIPE)
        q = len(q.stdout.read())
        print("")
        pbar = tqdm(desc="Jobs finished", total=int(runs), dynamic_ncols=True)
        while q > 2:
            t = subprocess.Popen("qstat", stdout=subprocess.PIPE)
            q = len(t.stdout.readlines())
            l = runs - (q - 2)
            pbar.n = int(l)
            pbar.refresh()
            time.sleep(2)
        else:
            print("\nDone processing, moving on to analysis")

    def shelx_write(self):
        shelxjob = open("shelxd_job.sh", "w")
        shelxjob.write("module load shelx\n")
        shelxjob.write("shelxd " + self.projname + "_fa")
        shelxjob.close()
        os.chmod("shelxd_job.sh", 0o775)

    def run_sagasu_proc(self):
        os.chdir(self.path)
        Path(self.projname).mkdir(parents=True, exist_ok=True)
        i = self.highres
        if self.clust == "l":
            tot = (self.lowres - self.highres) * (
                (self.highsites + 1) - self.lowsites
            )
            pbar = tqdm(desc="SHELXD", total=tot, dynamic_ncols=True)
        while not (i >= self.lowres):
            Path(os.path.join(self.projname, str(i))).mkdir(
                parents=True, exist_ok=True
            )
            i2 = i / 10
            j = self.highsites
            while not (j <= (self.lowsites - 1)):
                os.makedirs(
                    os.path.join(self.projname, str(i), str(j)), exist_ok=True
                )
                shutil.copy2(
                    self.insin, (os.path.join(self.projname, str(i), str(j)))
                )
                shutil.copy2(
                    self.hklin, (os.path.join(self.projname, str(i), str(j)))
                )
                shutil.copy2(
                    "shelxd_job.sh", (os.path.join(self.projname, str(i), str(j)))
                )
                workpath = os.path.join(self.path, self.projname, str(i), str(j))
                f = os.path.join(
                    self.path,
                    self.projname,
                    str(i),
                    str(j),
                    self.projname + "_fa.ins",
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
                    os.system(
                        "cd "
                        + workpath
                        + "; qsub -P i23 -q low.q -l mfree=8G -N sag_"
                        + str(i)
                        + "_"
                        + str(j)
                        + " -pe smp 40-33 -cwd ./shelxd_job.sh"
                    )
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
                print("Prepping for " + str(i2) + "Å, " + str(j) + " sites")
                csvfile = os.path.join(
                    self.path,
                    self.projname + "_results/" + str(i) + "_" + str(j) + ".csv",
                )
                numbers = str(i) + "_" + str(j)
                if self.clusteranalysis == "y":
                    print("***Bayesian Gaussian Mixture Analysis***")
                    #self.clustering_distance(csvfile, numbers)
                    clustering_distance_torun.append((csvfile, numbers))
                    print("***DBSCAN Analysis***")
                    #self.analysis(csvfile, numbers, i, j)
                    dbscan_torun.append((csvfile, numbers, i, j))
                    # print("***Generating Hexplots***")
                    # analysis_2(csvfile, numbers, i, j)
                    hexplots_torun.append((csvfile, numbers))
                else:
                    print("No cluster analysis requested")
                # print("***Outlier Analysis***")
                # ccalloutliers(csvfile, i, j)
                # ccweakoutliers(csvfile, i, j)
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
                # print("Results for " + str(i2) + "Å, " + str(j) + " sites:")
                csvfile = os.path.join(
                    self.path,
                    self.projname + "_results/" + str(i) + "_" + str(j) + ".csv",
                )
                numbers = str(i) + "_" + str(j)
                # print("***Generating ML Plots***")
                to_run_ML.append((csvfile, numbers))
                # plot_for_ML(csvfile, numbers)
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
        plt.scatter(df["CCWEAK"], df["CCALL"], c=df["PATFOM"], cmap='Blues', marker="o")
        #plt.axis("off")
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
        df.sort_values("CFOM", ascending=False, inplace=True, na_position='last')
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
        print("Best CFOM from", str(resolution/10), str(sitessearched), "is:", str(top_CFOM), end="\n")

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
        cmean = arr.mean(axis=0)
        csd = arr.std(axis=0)
        outliermask = ((arr[:, 0]) > (cmean[0] - (2 * csd[0]))) & (
            (arr[:, 1]) > (cmean[1] - (2 * csd[1]))
        )
        arr = arr[outliermask]
        mad = np.median(np.sqrt((arr[:, 0] - median) ** 2))
        ccallmax = heapq.nlargest(3, arr[:, 0])
        ccallmad = arr[:, 0] - median
        mad10 = sum(i > 10 * mad for i in ccallmad)
        mad9 = sum(i > 9 * mad for i in ccallmad)
        mad8 = sum(i > 8 * mad for i in ccallmad)
        mad7 = sum(i > 7 * mad for i in ccallmad)
        mad6 = sum(i > 6 * mad for i in ccallmad)
        mad5 = sum(i > 5 * mad for i in ccallmad)
        print("number of CCall with CCall - median > 10 * MAD = " + str(mad10))
        print("number of CCall with CCall - median >  9 * MAD = " + str(mad9))
        print("number of CCall with CCall - median >  8 * MAD = " + str(mad8))
        print("number of CCall with CCall - median >  7 * MAD = " + str(mad7))
        print("number of CCall with CCall - median >  6 * MAD = " + str(mad6))
        print("number of CCall with CCall - median >  5 * MAD = " + str(mad5))
        print("Three largest CCall values: " + str(ccallmax))
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
        ccweakmax = heapq.nlargest(3, arr[:, 1])
        ccweakmad = arr[:, 1] - median
        mad10 = sum(i > 10 * mad for i in ccweakmad)
        mad9 = sum(i > 9 * mad for i in ccweakmad)
        mad8 = sum(i > 8 * mad for i in ccweakmad)
        mad7 = sum(i > 7 * mad for i in ccweakmad)
        mad6 = sum(i > 6 * mad for i in ccweakmad)
        mad5 = sum(i > 5 * mad for i in ccweakmad)
        print("number of CCweak with CCweak - median > 10 * MAD = " + str(mad10))
        print("number of CCweak with CCweak - median >  9 * MAD = " + str(mad9))
        print("number of CCweak with CCweak - median >  8 * MAD = " + str(mad8))
        print("number of CCweak with CCweak - median >  7 * MAD = " + str(mad7))
        print("number of CCweak with CCweak - median >  6 * MAD = " + str(mad6))
        print("number of CCweak with CCweak - median >  5 * MAD = " + str(mad5))
        print("Three largest CCweak values: " + str(ccweakmax))
        print("")
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
        topall = str(top.reset_index(drop=True))
        print(
            """
        Here are the top 10 hits (ccall):

        """
        )
        print(topall)
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
        topweak = str(top.reset_index(drop=True))
        print(
            """
        Here are the top 10 hits (ccweak):

        """
        )
        print(topweak)
        cfom_df = pd.read_csv(self.projname + "_results/CFOM_PATFOM.csv", sep=",", names=["res", "sites", "CFOM", "PATFOM"])
        cfom_df.sort_values("CFOM", ascending=False, inplace=True)
        cfom_df["score"] = ((cfom_df["CFOM"]) - (((cfom_df["res"])) * (cfom_df["res"])) - (0.3 * (cfom_df["sites"])))
        top = cfom_df.head(10)
        top_CFOM = str(top.reset_index(drop=True))
        print(top_CFOM)
        with open('tophits.txt', 'w') as outfile:
            outfile.write(topall)
            outfile.write("\n")
            outfile.write(topweak)
            outfile.write("\n")
            outfile.write(top_CFOM)
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
            weak_df["res"], weak_df["sites"], weak_df["score"], cmap="viridis", edgecolor="none"
        )
        madplot = plt.gcf()
        madplot.savefig(self.projname + "_figures/ccweak.png", dpi=600)
        plt.clf()
        plt.cla()
        plt.close()
        ax = plt.axes(projection="3d")
        ax.plot_trisurf(cfom_df["res"], cfom_df["sites"], cfom_df["score"], cmap="viridis", edgecolor="none")
        madplot = plt.gcf()
        madplot.savefig(self.projname + "_figures/CFOM.png", dpi=600)
        # run phenix.emma on top 2 CFOM
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
                    sg = words[-1]
        print("The space group has been identified as " + sg)
        phenixemma = open("phenix_emma.sh", "w")
        phenixemma.write("module load phenix \n")
        phenixemma.write(
            "phenix.emma "
            + self.path
            + "/"
            + self.projname
            + "/"
            + str(firstres)
            + "/"
            + str(firstsites)
            + "/"
            + self.projname
            + "_fa.pdb "
            + self.path
            + "/"
            + self.projname
            + "/"
            + str(secondres)
            + "/"
            + str(secondsites)
            + "/"
            + self.projname
            + "_fa.pdb --tolerance=6 --space_group="
            + sg
            + " 2>&1 | tee -a phenixemma.log"
        )
        phenixemma.close()
        os.chmod("phenix_emma.sh", 0o775)
        os.system("./phenix_emma.sh")



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
            run.shelx_write(projname)
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
            #to_run = run.cleanup_prev()
            #pool.starmap(run.results, to_run)
            clustering_distance_torun, dbscan_torun, hexplots_torun, ccoutliers_torun = run.run_sagasu_analysis()
            #print("Clustering distance analysis...")
            #pool.starmap(run.clustering_distance, clustering_distance_torun)
            #print("DBScan")
            #pool.starmap(run.analysis, dbscan_torun)
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

