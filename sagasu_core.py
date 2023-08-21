#!/dls/science/groups/i23/pyenvs/sagasu_conda python3

from datetime import datetime
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import numpy as np
import pickle
import glob
import shutil
from pathlib import Path
import subprocess
import time
import requests
import paramiko
import json
from iotbx import mtz
from itertools import combinations




sns.set()


class core:
    def __init__(self):
        self.timestamp = datetime.now()
        self.path = os.getcwd()
        self.user = os.getenv('USER')

    def get_input(self):
        self.projname = input("Name of project: ")
        self.unitcell = str(input("Unit cell a b c al be ga: "))
        self.spacegroup = str(input("Spacegroup eg. P212121: "))
        # self.fa_path = input("Path to SHELXC outputs: ")
        self.fa_path = os.getcwd()
        self.highres = int(10 * float(input("High resolution cutoff for grid: ")))
        self.lowres = int(10 * float(input("Low resolution cutoff for grid: ")))
        self.highsites = int(input("Maximum number of sites to search: "))
        self.lowsites = int(input("Minimum number of sites to search: "))
        self.midsites = int(((self.highsites - self.lowsites) / 2) + self.lowsites)
        self.ntry = int(input("Number of trials: "))
        self.prasa_datain = input("HKL/mtz/sca input file for prasa: ")
        self.atomin = input("Anomalous scatterer: ")
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

    def get_slurm_token(self):
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        client.connect("wilson")
        command = "scontrol token lifespan=3600"
        stdin, stdout, stderr = client.exec_command(command)
        output = stdout.read().decode('utf-8')
        self.token = None
        for line in output.split('\n'):
            if line.startswith('SLURM_JWT'):
                self.token = line.split('=')[1].strip()
            else:
                pass
        client.close()

    def submit_shelxd_job_slurm(self):
        url = "http://slurm-rest.diamond.ac.uk:8443/slurm/v0.0.38/job/submit"
        headers = {'X-SLURM-USER-NAME': f'{self.user}', 'X-SLURM-USER-TOKEN': f'{self.token}', 'Content-Type': 'application/json'}


        script = "#!/bin/bash\n"
        script += "#SBATCH --job-name=sagasu\n"
        script += "#SBATCH --cpus-per-task=20\n"
        script += "#SBATCH --mem-per-cpu=1G\n"
        script += "#SBATCH --partition=cs04r\n"
        script += "#SBATCH --array=0-" + (str(len(self.shelxd_folder_paths) - 1)) + "\n"
        for i, folder in enumerate(self.shelxd_folder_paths):  
            script += f"#SBATCH --output={folder}/shelxd_output.log\n"
            script += f"#SBATCH --error={folder}/shelxd_error.log\n"
            script += f"/dls/science/groups/i23/scripts/chris/Sagasu/shelxd.sh {str(self.projname + '_fa')}\n"

        payload = {
            'script': script,
            'script_type': 'batch_script'
        }

        response = requests.post(url, headers=headers, data=json.dumps(payload))

        if response.status_code != 201:
            print("SHELXD: Failed to submit array job")
        else:
            print("SHELXD: Array job submitted")
            self.job_id_shelxd = response.json().get('job_id')


    def submit_afroprasa_job_slurm(self):
        url = "http://slurm-rest.diamond.ac.uk:8443/slurm/v0.0.38/job/submit"
        headers = {'X-SLURM-USER-NAME': f'{self.user}', 'X-SLURM-USER-TOKEN': f'{self.token}', 'Content-Type': 'application/json'}
        hr = str(self.highres / 10)
        lr = str(self.lowres / 10)

        script = "#!/bin/bash\n"
        script += "#SBATCH --job-name=sagasu\n"
        script += "#SBATCH --cpus-per-task=20\n"
        script += "#SBATCH --mem-per-cpu=1G\n"
        script += "#SBATCH --partition=cs04r\n"
        script += "#SBATCH --array=0-" + (str(len(self.prasa_folder_paths) - 1)) + "\n"
        for i, (folder, rescut) in enumerate(self.prasa_folder_paths):  
            script += f"#SBATCH --output={folder}/afroprasa_output.log\n"
            script += f"#SBATCH --error={folder}/afroprasa_error.log\n"
            script += f"/dls/science/groups/i23/scripts/chris/Sagasu/afroprasa.sh {self.atomin} {self.midsites} {str(rescut)} {self.ntry} {lr} {hr} {self.highsites} {self.lowsites}\n"

        payload = {
            'script': script,
            'script_type': 'batch_script'
        }

        response = requests.post(url, headers=headers, data=json.dumps(payload))

        if response.status_code != 201:
            print("Afroprasa: Failed to submit array job")
        else:
            print("Afroprasa: Array job submitted")     
            self.job_id_afroprasa = response.json().get('job_id')
  
       
    def wait_for_slurm_jobs(self):
        afroprasa_url = f"http://slurm-rest.diamond.ac.uk:8443/slurm/v0.0.38/jobs/{self.job_id_afroprasa}"
        shelxd_url = f"http://slurm-rest.diamond.ac.uk:8443/slurm/v0.0.38/jobs/{self.job_id_shelxd}" 
        headers = {'X-SLURM-USER-NAME': f'{self.user}', 'X-SLURM-USER-TOKEN': f'{self.token}'}
        afroprasa_finished = False
        shelxd_finished = False
        
        while not (afroprasa_finished and shelxd_finished):

            afroprasa_response = requests.get(afroprasa_url, headers=headers)
            shelxd_response = requests.get(shelxd_url, headers=headers)

            if afroprasa_response.status_code != 200 or shelxd_response.status_code != 200:
                print(f"Something went wrong. afroprasa status {afroprasa_response.text}, shelxd status {shelxd_response.text}")
                return False

            afroprasa_job_status = afroprasa_response.json().get('job_state')
            shelxd_job_status = shelxd_response.json().get('job_state')

            if not afroprasa_finished and afroprasa_job_status == 'COMPLETED':
                afroprasa_finished = True

            if not shelxd_finished and shelxd_job_status == 'COMPLETED':
                shelxd_finished = True

            if not (afroprasa_finished and shelxd_finished):
                time.sleep(5)

            if afroprasa_job_status in ['FAILED', 'CANCELLED'] or shelxd_job_status in ['FAILED', 'CANCELLED']:
                return False

        return True


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

    def shelxd_prep(self):
        os.system("module load ccp4")
        os.system(
            f"""
shelxc {self.projname} > /dev/null 2>&1 <<EOF
SAD aimless.sca
SFAC {self.atomin}
CELL {self.unitcell}
SPAG {self.spacegroup}
SHEL 999 {str(self.highres)}
FIND {str(self.lowsites)}
MIND -1.5
FRES 5
EOF
                  """
        )
# this preps for prasa and shelxd
    def prasa_prep(self):
        os.system("module load ccp4")
        if self.prasa_datain.endswith(".hkl" or ".HKL"):
            self.pointless = (
                "pointless HKLOUT pointless.mtz HKLIN "
                + str(self.prasa_datain)
                + " > /dev/null 2>&1"
            )
        elif self.prasa_datain.endswith(".sca" or ".SCA"):
            self.pointless = "pointless HKLOUT pointless.mtz SCAIN " + str(
                self.prasa_datain + " > /dev/null 2>&1"
            )
        else:
            print("\nMust be an mtz file in...\n")
            self.pointless = (
                f"pointless HKLOUT pointless.mtz HKLIN "
                + str(self.prasa_datain)
                + " > /dev/null 2>&1"
            )
        os.system(str(self.pointless))
        os.system(
            """
aimless HKLIN pointless.mtz HKLOUT aimless.mtz > /dev/null 2>&1 << eof
ANOMALOUS ON
eof
        """
        )
        os.system("mtz2sca aimless.mtz > /dev/null 2>&1")
        os.system(
            "ctruncate -hklin aimless.mtz -hklout truncate.mtz -colin '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]' > /dev/null 2>&1"
        )
        

    def get_unit_cell_and_sg(self):
        mtz = mtz.object("aimless.mtz")
        if mtz.crystals():
                crystal = mtz.crystals()[0]
                self.unitcell = crystal.unit_cell().parameters()
                self.spacegroup = crystal.space_group_info()
        else:
            pass
        if self.unitcell and self.spacegroup:
            pass
        else:
            self.spacegroup = input("Could not determine spacegroup, enter now (eg. P321): ")
            self.unitcell = input("Could not determine unit cell, enter now (eg. 150 150 45 90 90 120): ")

    def run_sagasu_proc(self):
        os.chdir(self.path)
        self.job_details = []
        self.shelxd_folder_paths = []
        self.prasa_folder_paths = []
        Path(self.projname).mkdir(parents=True, exist_ok=True)
        i = self.highres
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
                    os.chdir(self.path)
                elif self.clust == "c":    
                    self.shelxd_folder_paths.append(str(workpath))
                else:
                    print("error in input...")
                j = j - 1
            if self.clust == "c":
                os.makedirs(
                    os.path.join(self.projname, str(i), str(i) + "_prasa"),
                    exist_ok=True,
                )
                workpath = os.path.join(
                    self.path, self.projname, str(i), str(i) + "_prasa"
                )
                shutil.copy2(
                    "truncate.mtz",
                    (os.path.join(self.projname, str(i), str(i) + "_prasa")),
                )
                self.prasa_folder_paths.append((str(workpath), str(i/10)))
            else:
                pass
            i = i + 1

    def cleanup_prev(self):
        resultspath = os.path.join(self.path, self.projname + "_results")
        self.torun = []
        self.prasaruns = []
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
                    self.projname,
                    str(i),
                    str(j),
                    str(self.projname) + "_fa.lst",
                )
                self.torun.append((lstfile, i, j))
                j = j - 1
            prasaout = os.path.join(
                self.path, self.projname, str(i), str(i) + "_prasa", "prasa.txt"
            )
            self.prasaruns.append((prasaout, i))
            i = i + 1
        torun = self.torun
        prasaruns = self.prasaruns
        return torun, prasaruns

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

    def prasa_results(self, filename, i):
        # would be nice to find a way not to have two separate files...
        resfile = os.path.join(
            self.path, self.projname + "_results", "raw_prasa_" + str(i) + ".csv"
        )
        resfile_filt = os.path.join(
            self.path, self.projname + "_results", "prasa_" + str(i) + ".csv"
        )

        with open(filename, "r") as infile, open(resfile, "w") as outfile:
            filtered_file = ""
            for line in infile:
                if line.startswith("End of trial"):
                    if not line.endswith(" 0\n"):
                        outfile.write(line)
                else:
                    pass
        self.replace(resfile, "End of trial ", "")
        self.replace(resfile, "finalCC is ", "")
        self.replace(resfile, "CCrange is ", "")
        self.replace(resfile, "CCall is ", "")
        self.replace(resfile, "(candidate for a solution)", "")
        with open(resfile, "r") as file, open(resfile_filt, "w") as out:
            for line in file:
                lineout = line.replace(" ()", "")
                lineout = str(i) + ", " + lineout
                out.write(lineout)
        # might have to use the w.write(data[:-1]) to get rid of last whitespace in file

    def run_sagasu_analysis(self):
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
                ccoutliers_torun.append((csvfile, i, j))
                j = j - 1
            i = i + 1
        return ccoutliers_torun

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
        plt.scatter(df["CCWEAK"], df["CCALL"], marker="o")
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
        # outliermask = ((arr[:, 0]) > (cmean[0] - (2 * csd[0]))) & (
        #     (arr[:, 1]) > (cmean[1] - (2 * csd[1]))
        # )
        # chaging outlier mask to remove anything below the centroid in either direction
        outliermask = ((arr[:, 0]) > (cmean[0])) & ((arr[:, 1]) > (cmean[1]))
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
        # csd = arr.std(axis=0)
        # outliermask = ((arr[:, 0]) > (cmean[0] - (2 * csd[0]))) & (
        #     (arr[:, 1]) > (cmean[1] - (2 * csd[1]))
        # )
        outliermask = ((arr[:, 0]) > (cmean[0])) & ((arr[:, 1]) > (cmean[1]))
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

    def vectoroutliers_analysis(self, filename, resolution, sitessearched):
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
        pd.DataFrame.drop(
            df,
            labels=["linebeg", "CPUNO", "BEST", "TRY", "CFOM", "PATFOM"],
            axis=1,
            inplace=True,
        )
        ccall_mean = df["CCALL"].mean()
        ccweak_mean = df["CCWEAK"].mean()
        df["CCALL_VEC"] = df["CCALL"] - ccall_mean
        df["CCWEAK_VEC"] = df["CCWEAK"] - ccweak_mean
        df = df[df["CCALL_VEC"] > 0]
        df = df[df["CCWEAK_VEC"] > 0]
        df["VEC_DIFF"] = np.absolute(df["CCALL_VEC"] - df["CCWEAK_VEC"])
        df["COMB_VEC"] = np.sqrt(
            np.square(df["CCALL_VEC"]) + np.square(df["CCWEAK_VEC"])
        )
        df["NORM_VEC_DIFF"] = (df["VEC_DIFF"] / df["VEC_DIFF"].abs().max()) + 0.000001
        df["NORM_COMB_VEC"] = (df["COMB_VEC"] / df["COMB_VEC"].abs().max()) + 0.000001
        df["WEIGHTED"] = np.power(df["NORM_COMB_VEC"], 18) / np.cbrt(
            df["NORM_VEC_DIFF"]
        )
        df = df[df["WEIGHTED"] > 0.1]
        df["RES"] = resolution / 10
        df["SITES"] = sitessearched
        return df

    def vectoroutliers(self):
        all_data = pd.DataFrame()
        for resrange in range(self.highres, self.lowres):
            for siterange in range(self.lowsites, self.highsites):
                filename = os.path.join(
                    self.path,
                    self.projname
                    + "_results/"
                    + str(resrange)
                    + "_"
                    + str(siterange)
                    + ".csv",
                )
                data = self.vectoroutliers_analysis(filename, resrange, siterange)
                all_data = pd.concat([all_data, data], axis=0, ignore_index=True)
        all_data.sort_values(by=["COMB_VEC"], axis=0, inplace=True, ascending=False)
        customdata = np.stack(
            (all_data["RES"], all_data["SITES"], all_data["COMB_VEC"]), axis=1
        )
        fig = px.scatter(
            all_data,
            x="CCWEAK",
            y="CCALL",
            color="COMB_VEC",
            color_continuous_scale="Bluered_r",
        )
        hovertemplate = (
            "Res: %{customdata[0]} Å<br>"
            + "Sites: %{customdata[1]}<br>"
            + "Distance: %{customdata[2]:,.3f}<br>"
            + "CCWeak: %{x} <br>"
            + "CCAll: %{y}"
            + "<extra></extra>"
        )
        fig.update_traces(customdata=customdata, hovertemplate=hovertemplate)
        fig.write_html(self.projname + "_figures/vectoroutliers.html")

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
            "module load phenix && phenix.emma "
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

    def run_emma(self, emma_1, emma_2):
        command = f"module load phenix && phenix.emma --symmetry=scaled.mtz {str(file1)} {str(file2)}"  # Separate command and arguments
        result = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        return (emma_1, emma_2, result.stdout.decode('utf-8'))

    def get_filenames_for_emma(self):

        self.pdb_files_for_emma = [os.path.join(os.getcwd(), f"{self.projname}/{str(i)}/{str(j)}/{self.projname}_fa.pdb") for i in range(self.highres, self.lowres) for j in range(self.lowsites, self.highsites)] 
        parallel_filelist = list(combinations(pdb_files_for_emma, 2))
        return parallel_filelist

    def emma_correlation_plot(self, emma_results):
        pairs_pattern = r"Pairs:\s*(\d+)"
        singles_model1_pattern = r"Singles model 1:\s*(\d+)"
        singles_model2_pattern = r"Singles model 2:\s*(\d+)"        
        
        percentages = []
    
        for file in self.pdb_files_for_emma:
            filename = str(os.path.basename(os.path.dirname(os.path.dirname(file)))) + "_" + str(os.path.basename(os.path.dirname(file)))
            diagonalval = ([str(filename), str(filename), str(1)])
            percentages.append(diagonalval)

        for file1, file2, output in emma_results:
            pairs_match = re.search(pairs_pattern, output)
            singles_model1_match = re.search(singles_model1_pattern, output)
            singles_model2_match = re.search(singles_model2_pattern, output)

            if pairs_match and singles_model1_match and singles_model2_match:
                pairs_number = pairs_match.group(1)
                singles_model1_number = singles_model1_match.group(1)
                singles_model2_number = singles_model2_match.group(1)
            else:
                pass
            
            if pairs_number and singles_model1_number and singles_model2_number:
                percentage = np.around((float(pairs_number) / (float(pairs_number) + float(singles_model1_number) + float(singles_model2_number))), 2)
                filename1 = str(os.path.basename(os.path.dirname(os.path.dirname(file1)))) + "_" + str(os.path.basename(os.path.dirname(file1)))
                filename2 = str(os.path.basename(os.path.dirname(os.path.dirname(file2)))) + "_" + str(os.path.basename(os.path.dirname(file2)))
                percentages.append([str(filename1), str(filename2), str(percentage)])
            else:
                pass
        
        df = pd.DataFrame(percentages, columns=["filename_1", "filename_2", "percentage"])
        df_pivot1 = df.pivot(index='filename_2', columns='filename_1', values='percentage')
        df_pivot1.fillna(0, inplace=True)

        df_pivot2 = df.pivot(index='filename_1', columns='filename_2', values='percentage')
        df_pivot2.fillna(0, inplace=True)

        df_pivot3 = df_pivot1 + df_pivot2
        df_pivot3 = df_pivot3.multiply(100)
        df_pivot3.replace(200, 100, inplace=True)
        
        fig = px.imshow(df_pivot3, x=df_pivot3.columns, y=df_pivot3.index, text_auto=True, aspect="auto", color_continuous_scale="bluered")
        fig.update_layout(width=1000,height=1000)
        fig.write_html(self.projname + "_figures/emmamatrix.html")    
        
        
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
        
        <p><a href="./{projname}_figures/vectoroutliers.html">Vector Outliers Overview</a></p>
        
        <p><a href="./{projname}_figures/emmamatrix.html">Phenix EMMA Correlation Heatmap</a></p>
        
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
            for plot in sorted(
                glob.glob(os.path.join(self.path, (self.projname + "_figures"), "*ML*"))
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
