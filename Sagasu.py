# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 14:48:45 2020

@author: Chris
"""
#imports
#import subprocess
import os
from shutil import move
from tempfile import mkstemp
from os import remove

#environment setup
# =============================================================================
os.system("module load mx")
os.system("module load global/cluster")
# =============================================================================
path = os.getcwd()
print("You are here: " + path)
print("")

#user inputs
lowres = int((1 + (10 * float(input("Low resolution cutoff for grid: ")))))
highres = int((10 * float(input("High resolution cutoff for grid: "))))
lowsites = int(input("Minimum number of sites to search: "))
highsites = int(input("Maximum number of sites to search: "))
projname = ("wg99") #input("Name of project (SHELX prefix): ")
fa_path = ("/dls/i23/data/2020/nr23571-55/processing/wg99/xscale_9") #input("Path to SHELXC outputs: ")
clust = (str(input("Run on (c)luster or (l)ocal machine? ")))
ntry = int(input("Number of trials: "))
hklin = os.path.join(fa_path, projname + "_fa.hkl")
insin = os.path.join(fa_path, projname + "_fa.ins")

#write shelx job file
shelxjob = open("shelxd_job.sh","w")
shelxjob.write("module load shelx\n")
shelxjob.write("shelxd " + projname + "_fa")
shelxjob.close()
os.chmod('shelxd_job.sh', 0o775)

#replacement op
def replace(file_path, pattern, subst):
    file_path = os.path.abspath(file_path)
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    for line in old_file:
        new_file.write(subst if pattern in line else line)
    new_file.close()
    old_file.close()
    remove(file_path)
    move(abs_path, file_path)

#create subdirs,
os.system("mkdir "+path+"/"+projname)
for i in range(highres, lowres):
    os.system("mkdir "+projname+"/"+str(i))
    i2 = (i/10)
    for j in range(lowsites, highsites):
        os.system("mkdir "+projname+"/"+str(i)+"/"+str(j))
        os.system("cp "+insin+ " ./"+projname+"/"+str(i)+"/"+str(j))
        os.system("cp "+hklin+ " ./"+projname+"/"+str(i)+"/"+str(j))
        os.system("cp shelxd_job.sh "+" ./"+projname+"/"+str(i)+"/"+str(j))
        workpath = os.path.join(path, projname+"/"+str(i)+"/"+str(j))
        f = os.path.join(path, projname+"/"+str(i)+"/"+str(j)+"/"+projname+"_fa.ins")
        replace(f, "FIND", "FIND "+str(j)+"\n")
        replace(f, "SHEL", "SHEL 999 "+str(i2)+"\n")
        replace(f, "NTRY", "NTRY "+str(ntry)+"\n")
        if clust == 'l':
            print("Running on local machine, this may take some time...")
            os.system("cd "+workpath + "; ./shelxd_job.sh")
        elif clust == 'c':
            print("Submitting to the cluster, please run 'watch qstat' in another terminal to check on progress")
            os.system(" ")
    
print("Done. If nothing happened, make sure you pressed l or c at the cluster question.")