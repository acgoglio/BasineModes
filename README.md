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

