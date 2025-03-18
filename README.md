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

# while read LINE; do if [[ ${LINE:0:1} != '#' ]]; then echo $LINE ; MPython point_spt_diag.py $LINE ; fi ; done < idx_10pt.coo

