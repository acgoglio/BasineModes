# Prapare infiles
cdo selvar,MSL /data/inputs/METOCEAN/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/2014/01/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_msl

cdo setrtoc,0.9e+05,2e+05,1.01e+05 /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_msl /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_cost

cp /data/inputs/METOCEAN/historical/model/atmos/ECMWF/IFS_0125/analysis/6h/netcdf/2014/01/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/

cdo replace /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_cost /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_ok

mv -f /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc_ok /data/cmcc/ag15419/ECMWF_Med_modes/2014/12/20140101-ECMWF---AM0125-MEDATL-b20140102_an-fv10.00.nc 


# Compute Med Basin mode
#
# To run over an area:
# set and run run_basin_modes_cp.sh (it calls point_spt.py, which must be set too before running)
#
# If you just want to analyze a point with plots etc set and run point_spt_diag.py
# point_spt_diag2.py is the same of point_spt.py but includes diagnsotics
#
# NEW version:
f_point_area.py
run_basin_modes.py 

# SPython point_spt_diag_MSLP.py 200 80 prova
#
# while read LINE; do if [[ ${LINE:0:1} != '#' ]]; then echo $LINE ; MPython point_powspt_diag.py $LINE ; MPython point_ampspt_diag.py $LINE ; fi ; done < idx_10pt.coo


## To run on the whole domain
LPython run_basin_modes_amp_idx.py 300  421  0 128 1
LPython run_basin_modes_amp_idx.py 421  541  0 128 2
LPython run_basin_modes_amp_idx.py 541  661  0 128 3
LPython run_basin_modes_amp_idx.py 661  781  0 128 4
LPython run_basin_modes_amp_idx.py 781  901  0 128 5
LPython run_basin_modes_amp_idx.py 901  1021 0 128 6
LPython run_basin_modes_amp_idx.py 1021 1141 0 128 7
LPython run_basin_modes_amp_idx.py 1141 1261 0 128 8
LPython run_basin_modes_amp_idx.py 1261 1306 0 128 9
LPython run_basin_modes_amp_idx.py 300  421  128 255 10
LPython run_basin_modes_amp_idx.py 421  541  128 255 11
LPython run_basin_modes_amp_idx.py 541  661  128 255 12
LPython run_basin_modes_amp_idx.py 661  781  128 255 13
LPython run_basin_modes_amp_idx.py 781  901  128 255 14
LPython run_basin_modes_amp_idx.py 901  1021 128 255 15
LPython run_basin_modes_amp_idx.py 1021 1141 128 255 16
LPython run_basin_modes_amp_idx.py 1141 1261 128 255 17
LPython run_basin_modes_amp_idx.py 1261 1306 128 255 18
LPython run_basin_modes_amp_idx.py 300  421  255 379 19
LPython run_basin_modes_amp_idx.py 421  541  255 379 20
LPython run_basin_modes_amp_idx.py 541  661  255 379 21
LPython run_basin_modes_amp_idx.py 661  781  255 379 22
LPython run_basin_modes_amp_idx.py 781  901  255 379 23
LPython run_basin_modes_amp_idx.py 901  1021 255 379 24
LPython run_basin_modes_amp_idx.py 1021 1141 255 379 25
LPython run_basin_modes_amp_idx.py 1141 1261 255 379 26
LPython run_basin_modes_amp_idx.py 1261 1306 255 379 27

LPython run_basin_modes_pow_idx.py 300  421  0 128 1
LPython run_basin_modes_pow_idx.py 421  541  0 128 2
LPython run_basin_modes_pow_idx.py 541  661  0 128 3
LPython run_basin_modes_pow_idx.py 661  781  0 128 4
LPython run_basin_modes_pow_idx.py 781  901  0 128 5
LPython run_basin_modes_pow_idx.py 901  1021 0 128 6
LPython run_basin_modes_pow_idx.py 1021 1141 0 128 7
LPython run_basin_modes_pow_idx.py 1141 1261 0 128 8
LPython run_basin_modes_pow_idx.py 1261 1306 0 128 9
LPython run_basin_modes_pow_idx.py 300  421  128 255 10
LPython run_basin_modes_pow_idx.py 421  541  128 255 11
LPython run_basin_modes_pow_idx.py 541  661  128 255 12
LPython run_basin_modes_pow_idx.py 661  781  128 255 13
LPython run_basin_modes_pow_idx.py 781  901  128 255 14
LPython run_basin_modes_pow_idx.py 901  1021 128 255 15
LPython run_basin_modes_pow_idx.py 1021 1141 128 255 16
LPython run_basin_modes_pow_idx.py 1141 1261 128 255 17
LPython run_basin_modes_pow_idx.py 1261 1306 128 255 18
LPython run_basin_modes_pow_idx.py 300  421  255 379 19
LPython run_basin_modes_pow_idx.py 421  541  255 379 20
LPython run_basin_modes_pow_idx.py 541  661  255 379 21
LPython run_basin_modes_pow_idx.py 661  781  255 379 22
LPython run_basin_modes_pow_idx.py 781  901  255 379 23
LPython run_basin_modes_pow_idx.py 901  1021 255 379 24
LPython run_basin_modes_pow_idx.py 1021 1141 255 379 25
LPython run_basin_modes_pow_idx.py 1141 1261 255 379 26
LPython run_basin_modes_pow_idx.py 1261 1306 255 379 27
