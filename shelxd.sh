#!/bin/bash

# Script to  change resolution cutoff and number of sites for shelxd
#
# Aw 08/02/2016
# Updated Jul 2019 C Orr

# set project name
project="wzx"

# set path to processing directory
processpath=/dls/i23/data/2019/nr23571-10/processing/wzx/20190726/xscale/xscale_old_new/grid_wzx-old-new

# set resolution cutoffs (x10)
cutoff_high=30
cutoff_low=50

# set number of sites
sites_low=8
sites_high=22

############################################################################################################
# no changes below this line!!!
############################################################################################################

# set up shelxd script for qsub

echo module load shelx >> shelxd_job.sh
echo shelxd $project\_fa >> shelxd_job.sh
chmod +x shelxd_job.sh 

# load modules
module load shelx
module load global/cluster

if [ -d $processpath ]
  then
    echo
    echo "$processpath exists..."
  else
    echo "$processpath doesn't exist..."
    exit 0
fi

# 
# loop over different cutoffs 
#
cutoff=$cutoff_high 

while [ $cutoff -le $cutoff_low ]
do
  if [ -d $processpath/$project\_shelxd_$cutoff ]
  then
    echo "$processpath/$project\_shelxd_$cutoff exists..."
  else 
    mkdir $processpath/$project\_shelxd_$cutoff
  fi
  
  cd $processpath/$project\_shelxd_$cutoff
  
  # loop over number of sites
  sites=$sites_low 

  
    while [ $sites -le $sites_high ]
    do
      if [ -d $processpath/$project\_shelxd_$cutoff/sites_$sites ]
      then
        echo "$processpath/$project\_shelxd_$cutoff/sites_$sites exists..."
      else
        mkdir $processpath/$project\_shelxd_$cutoff/sites_$sites
      fi

      cd $processpath/$project\_shelxd_$cutoff/sites_$sites
      cp $processpath/scripts/$project\_fa.ins .
      cp $processpath/scripts/$project\_fa.hkl .
      res_cutoff=$(awk -v cutoff=$cutoff 'BEGIN { print cutoff / 10 }')
      sed -i "s/RESOLUTION/$res_cutoff/g" $project\_fa.ins
      sed -i "s/SITES/$sites/g" $project\_fa.ins
      echo `pwd`
      qsub -P i23 -q low.q -l h_vmem=8G -pe smp 20 -cwd $processpath/scripts/shelxd_job.sh  
      sites=$(( $sites + 1 ))
    done
  # end loop over number of sites

  cutoff=$(( $cutoff + 1 ))
done
# end loop over resolution
