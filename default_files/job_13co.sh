#!/bin/bash
#PBS -P iv23
#PBS -q express
#PBS -l walltime=14:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd
#PBS -N carver_13CO
#PBS -j oe
#PBS -m bea
#PBS -M u6645980@alumni.anu.edu.au

#to write levelpop_co.dat, include writepop keyword: radmc3d image npix 256 loadlambda fluxcons doppcatch inclline linelist nostar writepop sizepc 5.0 phi 0 incl 0 | tee output.txt
sizepc = 1.0
phi = 0
incl = 0
npix = 256

#get lines.inp file for 13CO
cp -p lines/lines_13co.inp .
mv lines_13co.inp lines.inp

#get temperature data for CO
cp -p temperatures/gas_temperature_CO.inp .
mv gas_temperature_CO.inp gas_temperature.inp

#CO Despotic
cp -p despotic/numberdens_13CO_despotic.inp .
mv numberdens_13CO_despotic.inp numberdens_13co.inp

#J=1-0
cp -p 13CO_J10/camera_wavelength_micron_13COJ10.inp .
mv camera_wavelength_micron_13COJ10.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_13CO_despotic_J10.txt
mv image.out image_13CO_despotic_J10.out
mv image_13CO_despotic_J10.out despotic/
mv output_13CO_despotic_J10.txt despotic/

#J=2-1
cp -p 13CO_J21/camera_wavelength_micron_13COJ21.inp .
mv camera_wavelength_micron_13COJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_13CO_despotic_J21.txt
mv image.out image_13CO_despotic_J21.out
mv image_13CO_despotic_J21.out despotic/
mv output_13CO_despotic_J21.txt despotic/

#J=3-2
cp -p 13CO_J32/camera_wavelength_micron_13COJ32.inp .
mv camera_wavelength_micron_13COJ32.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_13CO_despotic_J32.txt
mv image.out image_13CO_despotic_J32.out
mv image_13CO_despotic_J32.out despotic/
mv output_CO_despotic_J32.txt despotic/

#CO UCLCHEM
cp -p uclchem/numberdens_13CO_uclchem.inp .
mv numberdens_13CO_uclchem.inp numberdens_13co.inp

#J=1-0
cp -p 13CO_J10/camera_wavelength_micron_13COJ10.inp .
mv camera_wavelength_micron_13COJ10.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_13CO_uclchem_J10.txt
mv image.out image_13CO_uclchem_J10.out
mv image_13CO_uclchem_J10.out uclchem/
mv output_13CO_uclchem_J10.txt uclchem/

#J=2-1
cp -p 13CO_J21/camera_wavelength_micron_13COJ21.inp .
mv camera_wavelength_micron_13COJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_13CO_uclchem_J21.txt
mv image.out image_13CO_uclchem_J21.out
mv image_13CO_uclchem_J21.out uclchem/
mv output_13CO_uclchem_J21.txt uclchem/

#J=3-2
cp -p 13CO_J32/camera_wavelength_micron_13COJ32.inp .
mv camera_wavelength_micron_13COJ32.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_13CO_uclchem_J32.txt
mv image.out image_13CO_uclchem_J32.out
mv image_13CO_uclchem_J32.out uclchem/
mv output_13CO_uclchem_J32.txt uclchem/

rm radmc3d.out
rm camera_wavelength_micron*
rm numberdens_13co.inp
rm lines.inp
rm gas_temperature.inp
