# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 14:48:45 2020

@author: Chris
"""
#imports
import subprocess
import sys
import os

#environment setup
subprocess.run(["module load shelx"])
subprocess.run(["module load global/cluster"])
path = os.getcwd()
print("You are here: " + path)

#user inputs
lowres = input("Low resolution cutoff for grid: ")
highres = input("High resolution cutoff for grid: ")
lowsites = input("Minimum number of sites to search: ")
highsites = input("Maximum number of sites to search: ")
projname = input("Name of project (SHELX prefix): ")

