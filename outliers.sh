#/bin/bash

# Script to  change resolution cutoff and number of sites for shelxd
#
# Aw 08/02/2016
#

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


# check a few things first

if [ -d $processpath ]
  then
    echo
    echo "$processpath exists..."
  else
    echo "$processpath doesn't exist..."
    exit 0
fi

if [ -f $processpath/scripts/shelxd_outliers.py ]
  then 
    echo "$processpath/scripts/shelxd_outliers.py exists..."
  else
    echo "$processpath/scripts/shelxd_outliers.py doesn't exist..."
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
    echo "$processpath/$project\_shelxd_$cutoff doesn't exist..."
    exit 0
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
        echo "$processpath/$project\_shelxd_$cutoff/sites_$sites doesn't exist..."
        exit 0
      fi

      cd $processpath/$project\_shelxd_$cutoff/sites_$sites
      echo `pwd`
      $processpath/scripts/shelxd_outliers.py $project\_fa.lst > outliers.txt
      sites=$(( $sites + 1 ))
    done
  # end loop over number of sites

  cutoff=$(( $cutoff + 1 ))
done
# end loop over resolution
