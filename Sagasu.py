# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 14:48:45 2020

@author: Chris
"""
#imports
import subprocess
import os

#environment setup
# =============================================================================
# subprocess.run(["module load shelx"])
# subprocess.run(["module load global/cluster"])
# =============================================================================
path = os.getcwd()
print("You are here: " + path)
print("")

#user inputs
lowres = input("Low resolution cutoff for grid: ")
highres = input("High resolution cutoff for grid: ")
lowsites = input("Minimum number of sites to search: ")
highsites = input("Maximum number of sites to search: ")
projname = input("Name of project (SHELX prefix): ")
fa_path = input("Path to SHELXC outputs: ")
clust = input("Run on (c)luster or (l)ocal machine? ")
hklin = os.path.join(fa_path, projname + "_fa.hkl")
insin = os.path.join(fa_path, projname + "_fa.ins")

shelxjob = open("shelxd_job.sh","w")
shelxjob.write("module load shelx\n")
shelxjob.write("shelxd " + projname + "_fa")
shelxjob.close()
#subprocess.run(["chmod +x shelxd_job.sh"])

while (lowres >= highres):
    try:
        os.mkdir(path + "/" + projname + "/res_" + lowres)
        while (lowsites <= highsites):
            try:
                os.mkdir(path + "/" + projname + "/res_" + lowres + "/sites_" + lowsites)
                runfilepath = path + "/" + projname + "/res_" + lowres + "/sites_" + lowsites
                hklfilepath = os.path.join(runfilepath, projname + "_fa.hkl")
                origfafile = open(hklin, "r")
                fafile = open(hklfilepath, "w")
                fafile.write(origfafile)
                fafile.close()
                origfafile.close()
                lowsites = lowsites + 1
            except FileExistsError:
                pass
            lowres = lowres - 0.1
    except FileExistsError:
            pass                
else:
    print("Directory structure created...")