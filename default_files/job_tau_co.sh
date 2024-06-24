#!/bin/bash
#PBS -P jh2
#PBS -q normal
#PBS -l walltime=48:00:00
#PBS -l ncpus=4
#PBS -l mem=19GB
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd
#PBS -N carver_13CO_tau_5
#PBS -j oe
#PBS -m bea
#PBS -M u6645980@alumni.anu.edu.au

#to write levelpop_co.dat, include writepop keyword: radmc3d image npix 256 loadlambda fluxcons doppcatch inclline linelist nostar writepop sizepc 5.0 phi 0 incl 0 | tee output.txt
#REMEMBER sizepc should be 2*inputs.box_size you provide in inputs_gizmo_carver, since radmc wants it from left to right edge
sizepc=10.0
phi=0
incl=0
npix=256

#get lines.inp file for CO
cp -p lines/lines_co.inp .
mv lines_co.inp lines.inp

#get temperature data for CO
cp -p temperatures/gas_temperature_CO.inp .
mv gas_temperature_CO.inp gas_temperature.inp

#CO Despotic
cp -p despotic/numberdens_CO_despotic.inp .
mv numberdens_CO_despotic.inp numberdens_co.inp

#J=1-0
#cp -p CO_J10/camera_wavelength_micron_COJ10.inp .
#mv camera_wavelength_micron_COJ10.inp camera_wavelength_micron.inp
#radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl tracetau | tee tau_CO_despotic_J10.txt
#mv image.out tau_CO_despotic_J10.out
#mv tau_CO_despotic_J10.out despotic/
#mv tau_CO_despotic_J10.txt despotic/

#J=2-1
cp -p CO_J21/camera_wavelength_micron_COJ21.inp .
mv camera_wavelength_micron_COJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl tracetau | tee tau_CO_despotic_J21.txt
mv image.out tau_CO_despotic_J21.out
mv tau_CO_despotic_J21.out despotic/
mv tau_CO_despotic_J21.txt despotic/


#CO UCLCHEM
cp -p uclchem/numberdens_CO_uclchem.inp .
mv numberdens_CO_uclchem.inp numberdens_co.inp

#J=1-0
#cp -p CO_J10/camera_wavelength_micron_COJ10.inp .
#mv camera_wavelength_micron_COJ10.inp camera_wavelength_micron.inp
#radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl tracetau | tee tau_CO_uclchem_J10.txt
#mv image.out tau_CO_uclchem_J10.out
#mv tau_CO_uclchem_J10.out uclchem/
#mv tau_CO_uclchem_J10.txt uclchem/

#J=2-1
cp -p CO_J21/camera_wavelength_micron_COJ21.inp .
mv camera_wavelength_micron_COJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl tracetau | tee tau_CO_uclchem_J21.txt
mv image.out tau_CO_uclchem_J21.out
mv tau_CO_uclchem_J21.out uclchem/
mv tau_CO_uclchem_J21.txt uclchem/

rm radmc3d.out
rm camera_wavelength_micron*
rm numberdens_co.inp
rm lines.inp
rm gas_temperature.inp
