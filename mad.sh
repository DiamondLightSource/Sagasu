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


# a few checks
if [ -d $processpath ]
  then
    echo
    echo "$processpath exists..."
  else
    echo "$processpath doesn't exist..."
    exit 0
fi

if [ -d $processpath/results ]
  then
    echo
    echo "$processpath/results exists..."
  else
    mkdir $processpath/results
fi
# clean up first
if [ -f $processpath/results/mad5.txt ]
   then
     rm $processpath/results/mad5.txt
fi

if [ -f $processpath/results/mad6.txt ]
   then
     rm $processpath/results/mad6.txt
fi

if [ -f $processpath/results/mad7.txt ]
   then
     rm $processpath/results/mad7.txt
fi

if [ -f $processpath/results/mad8.txt ]
   then
     rm $processpath/results/mad8.txt
fi

if [ -f $processpath/results/mad9.txt ]
   then
     rm $processpath/results/mad9.txt
fi

if [ -f $processpath/results/mad10.txt ]
   then
     rm $processpath/results/mad10.txt
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
      awk '{if ($9=="5" && $12!="0")  {printf "cutoff: %3s sites: %3s hits: %3s \n", '"$cutoff"', '"$sites"', $12}}' outliers.txt >> $processpath/results/$project\_mad5.txt
      awk '{if ($9=="6" && $12!="0")  {printf "cutoff: %3s sites: %3s hits: %3s \n", '"$cutoff"', '"$sites"', $12}}' outliers.txt >> $processpath/results/$project\_mad6.txt
      awk '{if ($9=="7" && $12!="0")  {printf "cutoff: %3s sites: %3s hits: %3s \n", '"$cutoff"', '"$sites"', $12}}' outliers.txt >> $processpath/results/$project\_mad7.txt
      awk '{if ($9=="8" && $12!="0")  {printf "cutoff: %3s sites: %3s hits: %3s \n", '"$cutoff"', '"$sites"', $12}}' outliers.txt >> $processpath/results/$project\_mad8.txt
      awk '{if ($9=="9" && $12!="0")  {printf "cutoff: %3s sites: %3s hits: %3s \n", '"$cutoff"', '"$sites"', $12}}' outliers.txt >> $processpath/results/$project\_mad9.txt
      awk '{if ($9=="10" && $12!="0") {printf "cutoff: %3s sites: %3s hits: %3s \n", '"$cutoff"', '"$sites"', $12}}' outliers.txt >> $processpath/results/$project\_mad10.txt
      sites=$(( $sites + 1 ))
    done
  # end loop over number of sites

  cutoff=$(( $cutoff + 1 ))
done
# end loop over resolution
