#!/bin/bash
#PBS -P iv23
#PBS -q normal
#PBS -l walltime=14:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd
#PBS -N carver_C18O
#PBS -j oe
#PBS -m bea
#PBS -M u6645980@alumni.anu.edu.au

#to write levelpop_co.dat, include writepop keyword: radmc3d image npix 256 loadlambda fluxcons doppcatch inclline linelist nostar writepop sizepc 5.0 phi 0 incl 0 | tee output.txt
sizepc=5.0
phi=0
incl=0
npix=256

#get lines.inp file for C18O
cp -p lines/lines_c18o.inp .
mv lines_c18o.inp lines.inp

#get temperature data for CO
cp -p temperatures/gas_temperature_CO.inp .
mv gas_temperature_CO.inp gas_temperature.inp

#CO Despotic
cp -p despotic/numberdens_C18O_despotic.inp .
mv numberdens_C18O_despotic.inp numberdens_c18o.inp

#J=1-0
cp -p C18O_J10/camera_wavelength_micron_C18OJ10.inp .
mv camera_wavelength_micron_C18OJ10.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_C18O_despotic_J10.txt
mv image.out image_C18O_despotic_J10.out
mv image_C18O_despotic_J10.out despotic/
mv output_C18O_despotic_J10.txt despotic/

#J=2-1
cp -p C18O_J21/camera_wavelength_micron_C18OJ21.inp .
mv camera_wavelength_micron_C18OJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_C18O_despotic_J21.txt
mv image.out image_C18O_despotic_J21.out
mv image_C18O_despotic_J21.out despotic/
mv output_C18O_despotic_J21.txt despotic/

#J=3-2
cp -p C18O_J32/camera_wavelength_micron_C18OJ32.inp .
mv camera_wavelength_micron_C18OJ32.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_C18O_despotic_J32.txt
mv image.out image_C18O_despotic_J32.out
mv image_C18O_despotic_J32.out despotic/
mv output_C18O_despotic_J32.txt despotic/

#CO UCLCHEM
cp -p uclchem/numberdens_C18O_uclchem.inp .
mv numberdens_C18O_uclchem.inp numberdens_c18o.inp

#J=1-0
cp -p C18O_J10/camera_wavelength_micron_C18OJ10.inp .
mv camera_wavelength_micron_C18OJ10.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_C18O_uclchem_J10.txt
mv image.out image_C18O_uclchem_J10.out
mv image_C18O_uclchem_J10.out uclchem/
mv output_C18O_uclchem_J10.txt uclchem/

#J=2-1
cp -p C18O_J21/camera_wavelength_micron_C18OJ21.inp .
mv camera_wavelength_micron_C18OJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_C18O_uclchem_J21.txt
mv image.out image_C18O_uclchem_J21.out
mv image_C18O_uclchem_J21.out uclchem/
mv output_C18O_uclchem_J21.txt uclchem/

#J=3-2
cp -p C18O_J32/camera_wavelength_micron_C18OJ32.inp .
mv camera_wavelength_micron_C18OJ32.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_C18O_uclchem_J32.txt
mv image.out image_C18O_uclchem_J32.out
mv image_C18O_uclchem_J32.out uclchem/
mv output_C18O_uclchem_J32.txt uclchem/

rm radmc3d.out
rm camera_wavelength_micron*
rm numberdens_c18o.inp
rm lines.inp
rm gas_temperature.inp
