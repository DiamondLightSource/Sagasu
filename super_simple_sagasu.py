#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 16:53:17 2020

@author: vwg85559
"""
import os
import sagasu_core

pro_or_ana = input("p or a: ")
path = os.getcwd()
projname = input("proj: ")
fa_path = input("hkl2map path: ")
highres = input("highres: ")
lowres = input("lowres: ")
highsites = input("highsites: ")
lowsites = input("lowsites: ")
ntry = input("trys: ")
clust = 'c'

os.chdir(path)
clusteranalysis = "y"
projname = projname
pro_or_ana = str(pro_or_ana).lower()
highres = int((10 * float(highres)))
lowres = int((10 * float(lowres)))
highsites = int(highsites)
lowsites = int(lowsites)
insin = os.path.join(fa_path, projname + "_fa.ins")
hklin = os.path.join(fa_path, projname + "_fa.hkl")
ntry = int(ntry)
statusofrun = "-hold_jid "
clust = str(clust).lower()

if pro_or_ana == 'p':
    sagasu_core.writepickle(path,
                            projname,
                            lowres,
                            highres,
                            lowsites,
                            highsites,
                            ntry,
                            clusteranalysis
                            )
    sagasu_core.shelx_write(projname)
    sagasu_core.run_sagasu_proc(pro_or_ana,
                                projname,
                                highres,
                                lowres,
                                highsites,
                                lowsites,
                                insin,
                                hklin,
                                path,
                                ntry,
                                statusofrun,
                                clust
                                )

if pro_or_ana == 'a':
    sagasu_core.readpickle(path,
                           projname,
                           lowres,
                           highres,
                           lowsites,
                           highsites,
                           ntry,
                           clusteranalysis
                           )
    sagasu_core.cleanup_prev(path, projname, highres, lowres, highsites, lowsites)
    sagasu_core.for_ML_analysis(
        projname, highres, lowres, highsites, lowsites, path, clusteranalysis
        )