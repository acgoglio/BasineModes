#!/bin/bash
#
# by AC Goglio (CMCC)
# annachiara.goglio@cmcc.it
#
# Written: 06/02/2025
#
#set -u
set -e
#set -x
########################
# Inputs and Outputs

min_lat_idx=0
min_lon_idx=0
max_lat_idx=380
max_lon_idx=1307

work_dir='/work/cmcc/ag15419/basin_modes/'
infile_amppha='/work/cmcc/ag15419/basin_modes/basin_modes_ini.nc'
outfile='basin_modes.nc'

########################
echo "Load the env"
#    module load intel-2021.6.0/netcdf-c-threadsafe/4.9.0-25h5k
#    module load intel-2021.6.0/cdo-threadsafe/2.1.1-lyjsw
#    module load intel-2021.6.0/nco/5.0.6-jp6y4
#    module load intel-2021.6.0/ncview/2.1.8-sds5t
#    module load intel-2021.6.0/netcdf-c/4.9.0-cjqig
#    module load intel-2021.6.0/cdo-threadsafe/2.1.1-lyjsw
#    module load intel-2021.6.0/imagemagick/7.0.8-7-2475g
#    module load anaconda/3-2022.10
#    source activate /work/cmcc/ag15419/environment/mappyenv

########################
# Build the outfile 
cp -v ${infile_amppha} ${work_dir}/${outfile}

########################
echo "Loop on grid points:"
for lon_idx in $( seq $min_lon $max_lon ); do
    for lat_idx in $( seq $min_lat $max_lat ); do
           echo "Working on $lat_idx $lon_idx"
           #bsub -n 1 -q s_short -P 0510 -M 40G -o out -e err python point_spt.py $lat_idx $lon_idx ${work_dir}/${outfile}
    done          
    #sleep 10m
done
