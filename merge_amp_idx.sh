#!/bin/bash

# === Definizione intervalli di i (x) e j (y) ===
i_bounds=(300 420 540 660 780 900 1020 1140 1260 1307)
j_bounds=(0 127 254 380)

# === Step 1: creo file "vuoto" iniziale ===
echo "Initializing empty basin_modes_amp_med.nc..."
cp basin_modes_amp_1.nc basin_modes_amp_med.nc
cdo mulc,0 basin_modes_amp_med.nc basin_modes_amp_med.nc

count=1  # Contatore file

# === Step 2: ciclo su tutti i 27 box ===
for ((j=0; j<3; j++)); do
  for ((i=0; i<9; i++)); do
    echo "Processing basin_modes_amp_${count}.nc"

    xmin=${i_bounds[$i]}
    xmax=$(( ${i_bounds[$i+1]} - 1 ))
    ymin=${j_bounds[$j]}
    ymax=$(( ${j_bounds[$j+1]} - 1 ))

    # Estrai solo il box corretto (dall'originale)
    cdo selindexbox,$xmin,$xmax,$ymin,$ymax basin_modes_amp_${count}.nc box_${count}.nc

    # Rimuovi sossheig dal box estratto
    ncks -C -O -x -v sossheig box_${count}.nc box_${count}.nc

    # Somma il box ripulito al file finale
    cdo add basin_modes_amp_med.nc box_${count}.nc basin_modes_amp_med.nc

    # Elimina file temporaneo
    rm box_${count}.nc

    ((count++))
  done
done

# === Step 3: modifico attributi Phase -> Period e units deg -> hours ===
echo "Modifying standard_name and units in basin_modes_amp_med.nc..."

# Lista delle variabili Phase da modificare
phase_vars=(m0_T m1_T m2_T m3_T m4_T m5_T m6_T m7_T)

for var in "${phase_vars[@]}"; do
  ncatted -O \
    -a standard_name,"$var",m,c,"Period" \
    -a units,"$var",m,c,"hours" \
    basin_modes_amp_med.nc
done

echo "Operazione completata! File finale disponibile come basin_modes_amp_med.nc"

