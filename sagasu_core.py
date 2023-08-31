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
from iotbx.file_reader import any_file
from itertools import combinations


sns.set()


class core:
    def __init__(self):
        self.timestamp = datetime.now()
        self.path = os.getcwd()
        self.user = os.getenv("USER")
        self.env = {
            "CCP4I_TCLTK": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/bin",
            "CONDA_SHLVL": "2",
            "CONDA_EXE": "/dls_sw/apps/python/miniforge/4.10.0-0/bin/conda",
            "CLIBD_MON": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/lib/data/monomers/",
            "PYTHON_BASE_HOME": "/dls_sw/apps/python/miniforge/4.10.0-0",
            "LOADEDMODULES_modshare": "R/3.2.2:1:global/directories:1:shelx/ccp4:1:python/3.10:1:ccp4/8.0:1",
            "CONDA_PREFIX": "/dls/science/groups/i23/pyenvs/sagasu_conda",
            "MODULES_LMNOTUASKED_modshare": "R/3.2.2:1:global/directories:1:ccp4/8.0:1",
            "CCP4": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0",
            "CETC": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/etc",
            "CBIN": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/bin",
            "MODULES_CMD": "/usr/share/Modules/libexec/modulecmd.tcl",
            "MODULES_LMPREREQ": "ccp4/8.0&global/directories&R:shelx/ccp4&ccp4",
            "CONDA_PREFIX_1": "/dls_sw/apps/python/miniforge/4.10.0-0/envs/python3.10",
            "warpdoc": "/dls_sw/apps/ccp4/8.0.015/arp_warp_8.0/manual",
            "CLIB": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/lib",
            "CONDA_PYTHON_EXE": "/dls_sw/apps/python/miniforge/4.10.0-0/bin/python",
            "CRANK": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/share/ccp4i/crank",
            "CLIBD": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/lib/data",
            "CCP4_MASTER": "/dls_sw/apps/ccp4/8.0.015",
            "_LMFILES__modshare": "/dls_sw/apps/Modules/modulefiles/python/3.10:1:/dls_sw/apps/Modules/modulefiles/ccp4/8.0:1:/dls_sw/apps/Modules/modulefiles/R/3.2.2:1:/dls_sw/apps/Modules/modulefiles/global/directories:1:/dls_sw/apps/Modules/modulefiles/shelx/ccp4:1",
            "MODULES_LMALTNAME": "R/3.2.2&R/default&R:ccp4/8.0&ccp4/default&ccp4:shelx/ccp4&shelx/default&shelx:python/3.10&python/default&python",
            "GSETTINGS_SCHEMA_DIR": "/dls/science/groups/i23/pyenvs/sagasu_conda/share/glib-2.0/schemas",
            "MODULES_LMNOTUASKED": "global/directories:R/3.2.2:ccp4/8.0",
            "R_HOME": "/dls_sw/apps/R/3.2.2/lib64/R",
            "BALBES_ROOT": "/dls_sw/apps/ccp4/8.0.015/BALBES",
            "LOADEDMODULES": "global/directories:R/3.2.2:ccp4/8.0:shelx/ccp4:python/3.10",
            "MMCIFDIC": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/lib/ccp4/cif_mmdic.lib",
            "CHTML": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/html",
            "CONDA_PROMPT_MODIFIER": "(sagasu_conda)",
            "warpbin": "/dls_sw/apps/ccp4/8.0.015/arp_warp_8.0/bin/bin-x86_64-Linux",
            "SHELL": "/bin/bash",
            "CCP4I_TOP": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/share/ccp4i",
            "MANPATH_modshare": ":1:/var/cfengine/share/man:1",
            "MODULES_LMALTNAME_modshare": "R/3.2.2&R/default&R:1:shelx/ccp4&shelx/default&shelx:1:python/3.10&python/default&python:1:ccp4/8.0&ccp4/default&ccp4:1",
            "CINCL": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/include",
            "MANPATH": "/var/cfengine/share/man:",
            "MODULEPATH": "/etc/scl/modulefiles:/etc/scl/modulefiles:/etc/scl/modulefiles:/usr/share/Modules/modulefiles:/etc/modulefiles:/usr/share/modulefiles:/dls_sw/apps/Modules/modulefiles:/dls_sw/etc/modulefiles",
            "MODULES_LMPREREQ_modshare": "shelx/ccp4&ccp4:1:ccp4/8.0&global/directories&R:1",
            "MODULEPATH_modshare": "/etc/scl/modulefiles:1:/dls_sw/apps/Modules/modulefiles:1:/dls_sw/etc/modulefiles:1:/usr/share/Modules/modulefiles:2:/etc/modulefiles:2:/usr/share/modulefiles:2",
            "CEXAM": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/examples",
            "PATH": "/dls_sw/apps/python/miniforge/4.10.0-0/envs/python3.10/epics/bin/linux-x86_64:/dls/science/groups/i23/pyenvs/sagasu_conda/bin:/dls_sw/apps/python/miniforge/4.10.0-0/condabin:/dls_sw/apps/ccp4/8.0.015/arp_warp_8.0/bin/bin-x86_64-Linux:/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/etc:/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/bin:/dls_sw/apps/R/3.2.2/bin:/dls_sw/apps/R/3.2.2:/scratch/vwg85559/vscode/.vscode-server/bin/6c3e3dba23e8fadc360aed75ce363ba185c49794/bin/remote-cli:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/var/cfengine/bin:/home/i23user/bin:/home/vwg85559/bin/RIDL-master:/home/vwg85559/bin/SContent:/home/vwg85559/bin:/scratch/sw:/scratch/sw/pdbredo:/scratch/sw/bin:/scratch/sw/etc:/scratch/sw/usr/bin:/home/vwg85559/.local/bin:/home/i23user/bin/XZuiichi:/home/vwg85559/bin/cluster4xPrep:/scratch/Teams/usr/bin:/scratch/ChimeraX/usr/bin:/home/i23user/bin/Sagasu:/scratch/usr/lib64:/scratch/usr/share/doc:/home/i23user/bin/Sagasu:/scratch/blender-3.5.0-linux-x64/lib:/home/i23user/bin:/home/vwg85559/bin/RIDL-master:/home/vwg85559/bin/SContent:/home/vwg85559/bin:/scratch/sw:/scratch/sw/pdbredo:/scratch/sw/bin:/scratch/sw/etc:/scratch/sw/usr/bin:/home/vwg85559/.local/bin:/home/i23user/bin/XZuiichi:/home/vwg85559/bin/cluster4xPrep:/scratch/Teams/usr/bin:/scratch/ChimeraX/usr/bin:/home/i23user/bin/Sagasu:/scratch/usr/lib64:/scratch/usr/share/doc:/home/i23user/bin/Sagasu:/scratch/blender-3.5.0-linux-x64/lib:/home/i23user/bin:/home/vwg85559/bin/RIDL-master:/home/vwg85559/bin/SContent:/home/vwg85559/bin:/scratch/sw:/scratch/sw/pdbredo:/scratch/sw/bin:/scratch/sw/etc:/scratch/sw/usr/bin:/home/vwg85559/.local/bin:/home/i23user/bin/XZuiichi:/home/vwg85559/bin/cluster4xPrep:/scratch/Teams/usr/bin:/scratch/ChimeraX/usr/bin:/home/i23user/bin/Sagasu:/scratch/usr/lib64:/scratch/usr/share/doc:/home/i23user/bin/Sagasu:/scratch/blender-3.5.0-linux-x64/lib",
            "_LMFILES_": "/dls_sw/apps/Modules/modulefiles/global/directories:/dls_sw/apps/Modules/modulefiles/R/3.2.2:/dls_sw/apps/Modules/modulefiles/ccp4/8.0:/dls_sw/apps/Modules/modulefiles/shelx/ccp4:/dls_sw/apps/Modules/modulefiles/python/3.10",
            "MODULESHOME": "/usr/share/Modules",
            "CONDA_DEFAULT_ENV": "/dls/science/groups/i23/pyenvs/sagasu_conda",
        }

    def get_input(self):
        self.projname = input("Name of project: ")
        # self.unitcell = str(input("Unit cell a b c al be ga: "))
        # self.spacegroup = str(input("Spacegroup eg. P212121: "))
        # self.fa_path = input("Path to SHELXC outputs: ")
        self.prasa_datain = input("HKL/mtz/sca input file: ")
        self.get_unit_cell_and_sg()
        self.fa_path = os.getcwd()
        self.highres = int(10 * float(input("High resolution cutoff for grid: ")))
        self.lowres = int(10 * float(input("Low resolution cutoff for grid: ")))
        self.highsites = int(input("Maximum number of sites to search: "))
        self.lowsites = int(input("Minimum number of sites to search: "))
        self.midsites = int(((self.highsites - self.lowsites) / 2) + self.lowsites)
        self.ntry = int(input("Number of trials: "))
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
        output = stdout.read().decode("utf-8")
        self.token = None
        for line in output.split("\n"):
            if line.startswith("SLURM_JWT"):
                self.token = line.split("=")[1].strip()
            else:
                pass
        client.close()
        self.job_id_shelxd = []
        self.job_id_afroprasa = []

    def submit_shelxd_job_slurm(self):
        url = "https://slurm-rest.diamond.ac.uk:8443/slurm/v0.0.38/job/submit"
        headers = {
            "X-SLURM-USER-NAME": f"{self.user}",
            "X-SLURM-USER-TOKEN": self.token,
            "Content-Type": "application/json",
        }

        # if this doesnt work then use a tempfile workaround
        for i, folder in enumerate(self.shelxd_folder_paths):
            slurm_json = {
                "job": {
                    "name": "sagasu_shelxd",
                    "ntasks": 1,
                    "nodes": 1,
                    "current_working_directory": self.path,
                    "standard_input": "/dev/null",
                    "standard_output": f"{folder}/shelxd_output.log",
                    "standard_error": f"{folder}/shelxd_error.log",
                    "cpus_per_task": 20,
                    "partition": "cs04r",
                    "environment": self.env,
                },
                "script": f"#!/bin/bash\ncd {folder}\n/dls/science/groups/i23/scripts/chris/Sagasu_slurm/Sagasu/shelxd.sh {str(self.projname + '_fa')}",
            }

            response = requests.post(url, headers=headers, data=json.dumps(slurm_json))

            if response.status_code != 200:
                print(response.text)
            else:
                self.job_id_shelxd.append(response.json().get("job_id"))

    def submit_afroprasa_job_slurm(self):
        url = "https://slurm-rest.diamond.ac.uk:8443/slurm/v0.0.38/job/submit"
        headers = {
            "X-SLURM-USER-NAME": f"{self.user}",
            "X-SLURM-USER-TOKEN": f"{self.token}",
            "Content-Type": "application/json",
        }
        hr = str(self.highres / 10)
        lr = str(self.lowres / 10)

        # if this doesnt work then use a tempfile workaround
        for i, (folder, rescut) in enumerate(self.prasa_folder_paths):
            slurm_json = {
                "job": {
                    "name": "sagasu_afroprasa",
                    "ntasks": 1,
                    "nodes": 1,
                    "current_working_directory": self.path,
                    "standard_input": "/dev/null",
                    "standard_output": f"{folder}/afroprasa_output.log",
                    "standard_error": f"{folder}/afroprasa_error.log",
                    "cpus_per_task": 20,
                    "partition": "cs04r",
                    "environment": self.env,
                },
                "script": f"#!/bin/bash\ncd {folder}\n/dls/science/groups/i23/scripts/chris/Sagasu_slurm/Sagasu/afroprasa.sh {self.atomin} {self.midsites} {str(rescut)} {self.ntry} {lr} {hr} {self.highsites} {self.lowsites}",
            }

            response = requests.post(url, headers=headers, data=json.dumps(slurm_json))

            if response.status_code == 200:
                self.job_id_afroprasa.append(response.json().get("job_id"))
            else:
                print(response.status_code)
                print(response.text)

    def wait_for_slurm_jobs(self):
        jobs_url = "https://slurm-rest.diamond.ac.uk:8443/slurm/v0.0.38/job/"
        headers = {
            "X-SLURM-USER-NAME": f"{self.user}",
            "X-SLURM-USER-TOKEN": f"{self.token}",
        }

        shelxd_finished = False
        shelxd_unfinished = []
        afroprasa_finished = False
        afroprasa_unfinished = []

        while not (afroprasa_finished and shelxd_finished):
            if not shelxd_finished:
                for shelxd_job in self.job_id_shelxd:
                    shelxd_response = requests.get(
                        f"{jobs_url}/{shelxd_job}", headers=headers
                    )
                    if shelxd_response.status_code == 200:
                        if not shelxd_response.json().get("job_state") == "COMPLETED":
                            shelxd_unfinished.append(shelxd_job)
                        elif (
                            shelxd_response.json().get("job_state") == "FAILED"
                            or shelxd_response.json().get("job_state") == "CANCELLED"
                        ):
                            print(
                                "Something went wrong and the job failed or was cancelled"
                            )
                    else:
                        print(shelxd_response.status_code)
                        print(shelxd_response.text)
                    self.job_id_shelxd = shelxd_unfinished

            if not afroprasa_finished:
                for afroprasa_job in self.job_id_afroprasa:
                    afroprasa_response = requests.get(
                        f"{jobs_url}/{afroprasa_job}", headers=headers
                    )
                    if afroprasa_response.status_code == 200:
                        if (
                            not afroprasa_response.json().get("job_state")
                            == "COMPLETED"
                        ):
                            afroprasa_unfinished.append(afroprasa_job)
                        elif (
                            afroprasa_response.json().get("job_state") == "FAILED"
                            or afroprasa_response.json().get("job_state") == "CANCELLED"
                        ):
                            print(
                                "Something went wrong and the job failed or was cancelled"
                            )
                    else:
                        print(afroprasa_response.status_code)
                        print(afroprasa_response.text)
                self.job_id_shelxd = afroprasa_unfinished

            if len(self.job_id_shelxd) == 0:
                shelxd_finished = True
            if len(self.job_id_afroprasa) == 0:
                afroprasa_finished = True

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
                    self.atomin,
                    self.prasa_datain,
                    self.midsites,
                    self.unitcell,
                    self.spacegroup,
                ],
                f,
            )

    def readpickle(self):
        with open("inps.pkl", "rb") as f:
            pickle.load(f)
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
                self.atomin,
                self.prasa_datain,
                self.midsites,
                self.unitcell,
                self.spacegroup,
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
        print("Running SHEXC...")
        os.system("module load ccp4")
        print("Loaded shelx")
        print(
            f"""
shelxc {self.projname} > /dev/null 2>&1 <<EOF
SAD aimless.sca
SFAC {(self.atomin).upper()}
CELL {self.unitcell}
SPAG {self.spacegroup}
SHEL 999 {str(self.highres)}
FIND {str(self.lowsites)}
MIND -1.5
FRES 5
EOF
                  """
        )
        os.system(
            f"""
shelxc {self.projname} > /dev/null 2>&1 <<EOF
SAD aimless.sca
SFAC {(self.atomin).upper()}
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
        try:
            read_data_file = any_file(self.prasa_datain)
            data_file_symm = read_data_file.crystal_symmetry()
            symm_as_py_code = data_file_symm.as_py_code()
            unit_cell_match = re.search(r"unit_cell=\((.*?)\)", symm_as_py_code)
            self.unitcell = unit_cell_match.group(1)
            self.unitcell = self.unitcell.replace(",", "")
            space_group_match = re.search(
                r'space_group_symbol="([^"]+)"', symm_as_py_code
            )
            self.spacegroup = space_group_match.group(1)
            self.spacegroup = self.spacegroup.replace(" ", "")
        except:
            pass
        if self.unitcell and self.spacegroup:
            print(
                f"Spacegroup and unit cell identified as {str(self.spacegroup)}, {str(self.unitcell)}"
            )
        else:
            self.spacegroup = input(
                "Could not determine spacegroup, enter now (eg. P321): "
            )
            self.unitcell = input(
                "Could not determine unit cell, enter now (eg. 150 150 45 90 90 120): "
            )

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
                self.prasa_folder_paths.append((str(workpath), str(i / 10)))
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
        return (emma_1, emma_2, result.stdout.decode("utf-8"))

    def get_filenames_for_emma(self):
        self.pdb_files_for_emma = [
            os.path.join(
                os.getcwd(), f"{self.projname}/{str(i)}/{str(j)}/{self.projname}_fa.pdb"
            )
            for i in range(self.highres, self.lowres)
            for j in range(self.lowsites, self.highsites)
        ]
        parallel_filelist = list(combinations(self.pdb_files_for_emma, 2))
        return parallel_filelist

    def emma_correlation_plot(self, emma_results):
        pairs_pattern = r"Pairs:\s*(\d+)"
        singles_model1_pattern = r"Singles model 1:\s*(\d+)"
        singles_model2_pattern = r"Singles model 2:\s*(\d+)"

        percentages = []

        for file in self.pdb_files_for_emma:
            filename = (
                str(os.path.basename(os.path.dirname(os.path.dirname(file))))
                + "_"
                + str(os.path.basename(os.path.dirname(file)))
            )
            diagonalval = [str(filename), str(filename), str(1)]
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
                percentage = np.around(
                    (
                        float(pairs_number)
                        / (
                            float(pairs_number)
                            + float(singles_model1_number)
                            + float(singles_model2_number)
                        )
                    ),
                    2,
                )
                filename1 = (
                    str(os.path.basename(os.path.dirname(os.path.dirname(file1))))
                    + "_"
                    + str(os.path.basename(os.path.dirname(file1)))
                )
                filename2 = (
                    str(os.path.basename(os.path.dirname(os.path.dirname(file2))))
                    + "_"
                    + str(os.path.basename(os.path.dirname(file2)))
                )
                percentages.append([str(filename1), str(filename2), str(percentage)])
            else:
                pass

        df = pd.DataFrame(
            percentages, columns=["filename_1", "filename_2", "percentage"]
        )
        df_pivot1 = df.pivot(
            index="filename_2", columns="filename_1", values="percentage"
        )
        df_pivot1.fillna(0, inplace=True)

        df_pivot2 = df.pivot(
            index="filename_1", columns="filename_2", values="percentage"
        )
        df_pivot2.fillna(0, inplace=True)

        df_pivot3 = df_pivot1 + df_pivot2
        df_pivot3 = df_pivot3.multiply(100)
        df_pivot3.replace(200, 100, inplace=True)

        fig = px.imshow(
            df_pivot3,
            x=df_pivot3.columns,
            y=df_pivot3.index,
            text_auto=True,
            aspect="auto",
            color_continuous_scale="bluered",
        )
        fig.update_layout(width=1000, height=1000)
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


env = {
    "CCP4I_TCLTK": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/bin",
    "CONDA_SHLVL": "2",
    "CONDA_EXE": "/dls_sw/apps/python/miniforge/4.10.0-0/bin/conda",
    "CLIBD_MON": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/lib/data/monomers/",
    "PYTHON_BASE_HOME": "/dls_sw/apps/python/miniforge/4.10.0-0",
    "MODULES_RUN_QUARANTINE": "LD_LIBRARY_PATH LD_PRELOAD",
    "LANG": "en_GB.UTF-8",
    "HISTCONTROL": "ignoredups",
    "HOSTNAME": "ws370.diamond.ac.uk",
    "OPENBLAS_NUM_THREADS": "32",
    "softwaredir": "/dls_sw/apps",
    "PATH_modshare": "/scratch/usr/lib64:1:/scratch/ChimeraX/usr/bin:1:/scratch/sw:1:/dls_sw/apps/ccp4/8.0.015/arp_warp_8.0/bin/bin-x86_64-Linux:1:/scratch/sw/usr/bin:1:/dls_sw/apps/R/3.2.2:1:/usr/bin:1:/home/vwg85559/.local/bin:1:/scratch/sw/etc:1:/home/vwg85559/bin:1:/usr/share/Modules/bin:1:/usr/local/bin:1:/scratch/vwg85559/vscode/.vscode-server/bin/6c3e3dba23e8fadc360aed75ce363ba185c49794/bin/remote-cli:1:/var/cfengine/bin:1:/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/etc:1:/scratch/sw/pdbredo:1:/home/vwg85559/bin/cluster4xPrep:1:/home/i23user/bin/Sagasu:1:/home/vwg85559/bin/SContent:1:/scratch/sw/bin:1:/home/i23user/bin/XZuiichi:1:/dls_sw/apps/R/3.2.2/bin:1:/scratch/Teams/usr/bin:1:/home/i23user/bin:1:/usr/sbin:1:/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/bin:1:/scratch/blender-3.5.0-linux-x64/lib:1:/usr/local/sbin:1:/home/vwg85559/bin/RIDL-master:1:/scratch/usr/share/doc:1",
    "LOADEDMODULES_modshare": "R/3.2.2:1:global/directories:1:shelx/ccp4:1:python/3.10:1:ccp4/8.0:1",
    "COLORTERM": "truecolor",
    "GSETTINGS_SCHEMA_DIR_CONDA_BACKUP": "",
    "CONDA_PREFIX": "/dls/science/groups/i23/pyenvs/sagasu_conda",
    "MODULES_LMNOTUASKED_modshare": "R/3.2.2:1:global/directories:1:ccp4/8.0:1",
    "ISPYB_CONFIG_FILE": "/dls_sw/dasc/mariadb/credentials/ispyb.cfg",
    "VSCODE_GIT_ASKPASS_EXTRA_ARGS": "",
    "CCP4": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0",
    "CCP4_OPEN": "UNKNOWN",
    "CETC": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/etc",
    "S_COLORS": "auto",
    "CCP4_SCR": "/dls/tmp/vwg85559",
    "_CE_M": "",
    "which_declare": "declare -f",
    "CBIN": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/bin",
    "XDG_SESSION_ID": "65",
    "MODULES_CMD": "/usr/share/Modules/libexec/modulecmd.tcl",
    "USER": "vwg85559",
    "MODULES_LMPREREQ": "ccp4/8.0&global/directories&R:shelx/ccp4&ccp4",
    "CONDA_PREFIX_1": "/dls_sw/apps/python/miniforge/4.10.0-0/envs/python3.10",
    "warpdoc": "/dls_sw/apps/ccp4/8.0.015/arp_warp_8.0/manual",
    "SELINUX_ROLE_REQUESTED": "",
    "GFORTRAN_UNBUFFERED_PRECONNECTED": "1",
    "CLIB": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/lib",
    "PWD": "/dls/i23/data/2023/cm33851-3/processing/chris/slurmy/20/20_prasa",
    "SSH_ASKPASS": "/usr/libexec/openssh/gnome-ssh-askpass",
    "HOME": "/home/vwg85559",
    "CONDA_PYTHON_EXE": "/dls_sw/apps/python/miniforge/4.10.0-0/bin/python",
    "BROWSER": "/scratch/vwg85559/vscode/.vscode-server/bin/6c3e3dba23e8fadc360aed75ce363ba185c49794/bin/helpers/browser.sh",
    "VSCODE_GIT_ASKPASS_NODE": "/scratch/vwg85559/vscode/.vscode-server/bin/6c3e3dba23e8fadc360aed75ce363ba185c49794/node",
    "TERM_PROGRAM": "vscode",
    "SSH_CLIENT": "172.23.100.100 36600 22",
    "TERM_PROGRAM_VERSION": "1.81.1",
    "SELINUX_LEVEL_REQUESTED": "",
    "CRANK": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/share/ccp4i/crank",
    "CLIBD": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/lib/data",
    "XDG_DATA_DIRS": "/home/vwg85559/.local/share/flatpak/exports/share:/var/lib/flatpak/exports/share:/usr/local/share:/usr/share",
    "CCP4_MASTER": "/dls_sw/apps/ccp4/8.0.015",
    "_LMFILES__modshare": "/dls_sw/apps/Modules/modulefiles/python/3.10:1:/dls_sw/apps/Modules/modulefiles/ccp4/8.0:1:/dls_sw/apps/Modules/modulefiles/R/3.2.2:1:/dls_sw/apps/Modules/modulefiles/global/directories:1:/dls_sw/apps/Modules/modulefiles/shelx/ccp4:1",
    "MODULES_LMALTNAME": "R/3.2.2&R/default&R:ccp4/8.0&ccp4/default&ccp4:shelx/ccp4&shelx/default&shelx:python/3.10&python/default&python",
    "_CE_CONDA": "",
    "GSETTINGS_SCHEMA_DIR": "/dls/science/groups/i23/pyenvs/sagasu_conda/share/glib-2.0/schemas",
    "MODULES_LMNOTUASKED": "global/directories:R/3.2.2:ccp4/8.0",
    "VSCODE_IPC_HOOK_CLI": "/run/user/1015129/vscode-ipc-dc87649e-146a-45f3-bc9a-daf0d0a9465a.sock",
    "R_HOME": "/dls_sw/apps/R/3.2.2/lib64/R",
    "BALBES_ROOT": "/dls_sw/apps/ccp4/8.0.015/BALBES",
    "LOADEDMODULES": "global/directories:R/3.2.2:ccp4/8.0:shelx/ccp4:python/3.10",
    "MMCIFDIC": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/lib/ccp4/cif_mmdic.lib",
    "CHTML": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/html",
    "CONDA_PROMPT_MODIFIER": "(sagasu_conda)",
    "warpbin": "/dls_sw/apps/ccp4/8.0.015/arp_warp_8.0/bin/bin-x86_64-Linux",
    "MAIL": "/var/spool/mail/vwg85559",
    "VSCODE_GIT_ASKPASS_MAIN": "/scratch/vwg85559/vscode/.vscode-server/bin/6c3e3dba23e8fadc360aed75ce363ba185c49794/extensions/git/dist/askpass-main.js",
    "localTEMP": "/tmp/vwg85559",
    "networkTEMP": "/dls/tmp/vwg85559",
    "TERM": "xterm-256color",
    "SHELL": "/bin/bash",
    "CCP4I_TOP": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/share/ccp4i",
    "MANPATH_modshare": ":1:/var/cfengine/share/man:1",
    "MODULES_LMALTNAME_modshare": "R/3.2.2&R/default&R:1:shelx/ccp4&shelx/default&shelx:1:python/3.10&python/default&python:1:ccp4/8.0&ccp4/default&ccp4:1",
    "CINCL": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/include",
    "SELINUX_USE_CURRENT_RANGE": "",
    "SHLVL": "5",
    "VSCODE_GIT_IPC_HANDLE": "/run/user/1015129/vscode-git-22e00c5660.sock",
    "MANPATH": "/var/cfengine/share/man:",
    "MODULEPATH": "/etc/scl/modulefiles:/etc/scl/modulefiles:/etc/scl/modulefiles:/usr/share/Modules/modulefiles:/etc/modulefiles:/usr/share/modulefiles:/dls_sw/apps/Modules/modulefiles:/dls_sw/etc/modulefiles",
    "LOGNAME": "vwg85559",
    "MODULES_LMPREREQ_modshare": "shelx/ccp4&ccp4:1:ccp4/8.0&global/directories&R:1",
    "DBUS_SESSION_BUS_ADDRESS": "unix:path=/run/user/1015129/bus",
    "GIT_ASKPASS": "/scratch/vwg85559/vscode/.vscode-server/bin/6c3e3dba23e8fadc360aed75ce363ba185c49794/extensions/git/dist/askpass.sh",
    "MODULEPATH_modshare": "/etc/scl/modulefiles:1:/dls_sw/apps/Modules/modulefiles:1:/dls_sw/etc/modulefiles:1:/usr/share/Modules/modulefiles:2:/etc/modulefiles:2:/usr/share/modulefiles:2",
    "CEXAM": "/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/examples",
    "PATH": "/dls_sw/apps/python/miniforge/4.10.0-0/envs/python3.10/epics/bin/linux-x86_64:/dls/science/groups/i23/pyenvs/sagasu_conda/bin:/dls_sw/apps/python/miniforge/4.10.0-0/condabin:/dls_sw/apps/ccp4/8.0.015/arp_warp_8.0/bin/bin-x86_64-Linux:/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/etc:/dls_sw/apps/ccp4/8.0.015/ccp4-8.0/bin:/dls_sw/apps/R/3.2.2/bin:/dls_sw/apps/R/3.2.2:/scratch/vwg85559/vscode/.vscode-server/bin/6c3e3dba23e8fadc360aed75ce363ba185c49794/bin/remote-cli:/usr/share/Modules/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/var/cfengine/bin:/home/i23user/bin:/home/vwg85559/bin/RIDL-master:/home/vwg85559/bin/SContent:/home/vwg85559/bin:/scratch/sw:/scratch/sw/pdbredo:/scratch/sw/bin:/scratch/sw/etc:/scratch/sw/usr/bin:/home/vwg85559/.local/bin:/home/i23user/bin/XZuiichi:/home/vwg85559/bin/cluster4xPrep:/scratch/Teams/usr/bin:/scratch/ChimeraX/usr/bin:/home/i23user/bin/Sagasu:/scratch/usr/lib64:/scratch/usr/share/doc:/home/i23user/bin/Sagasu:/scratch/blender-3.5.0-linux-x64/lib:/home/i23user/bin:/home/vwg85559/bin/RIDL-master:/home/vwg85559/bin/SContent:/home/vwg85559/bin:/scratch/sw:/scratch/sw/pdbredo:/scratch/sw/bin:/scratch/sw/etc:/scratch/sw/usr/bin:/home/vwg85559/.local/bin:/home/i23user/bin/XZuiichi:/home/vwg85559/bin/cluster4xPrep:/scratch/Teams/usr/bin:/scratch/ChimeraX/usr/bin:/home/i23user/bin/Sagasu:/scratch/usr/lib64:/scratch/usr/share/doc:/home/i23user/bin/Sagasu:/scratch/blender-3.5.0-linux-x64/lib:/home/i23user/bin:/home/vwg85559/bin/RIDL-master:/home/vwg85559/bin/SContent:/home/vwg85559/bin:/scratch/sw:/scratch/sw/pdbredo:/scratch/sw/bin:/scratch/sw/etc:/scratch/sw/usr/bin:/home/vwg85559/.local/bin:/home/i23user/bin/XZuiichi:/home/vwg85559/bin/cluster4xPrep:/scratch/Teams/usr/bin:/scratch/ChimeraX/usr/bin:/home/i23user/bin/Sagasu:/scratch/usr/lib64:/scratch/usr/share/doc:/home/i23user/bin/Sagasu:/scratch/blender-3.5.0-linux-x64/lib",
    "_LMFILES_": "/dls_sw/apps/Modules/modulefiles/global/directories:/dls_sw/apps/Modules/modulefiles/R/3.2.2:/dls_sw/apps/Modules/modulefiles/ccp4/8.0:/dls_sw/apps/Modules/modulefiles/shelx/ccp4:/dls_sw/apps/Modules/modulefiles/python/3.10",
    "MODULESHOME": "/usr/share/Modules",
    "CONDA_DEFAULT_ENV": "/dls/science/groups/i23/pyenvs/sagasu_conda",
    "HISTSIZE": "1000",
    "XML_CATALOG_FILES": "file:///dls/science/groups/i23/pyenvs/sagasu_conda/etc/xml/catalog file:///etc/xml/catalog",
}
