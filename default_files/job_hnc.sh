#!/bin/bash
#PBS -P iv23
#PBS -q express
#PBS -l walltime=18:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd
#PBS -N carver_HNC
#PBS -j oe
#PBS -m bea
#PBS -M u6645980@alumni.anu.edu.au

#to write levelpop_co.dat, include writepop keyword: radmc3d image npix 256 loadlambda fluxcons doppcatch inclline linelist nostar writepop sizepc 5.0 phi 0 incl 0 | tee output.txt
sizepc=5.0
phi=0
incl=0
npix=256

#get lines.inp file for HNC
cp -p lines/lines_hnc.inp .
mv lines_hnc.inp lines.inp

#get temperature data for HNC
cp -p temperatures/gas_temperature_HNC.inp .
mv gas_temperature_HNC.inp gas_temperature.inp

#HNC UCLCHEM
cp -p uclchem/numberdens_HNC_uclchem.inp .
mv numberdens_HNC_uclchem.inp numberdens_hnc.inp

#J=1-0
cp -p HNC_J10/camera_wavelength_micron_HNCJ10.inp .
mv camera_wavelength_micron_HNCJ10.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl setthreads $PBS_NCPUS | tee output_HNC_uclchem_J10.txt
mv image.out image_HNC_uclchem_J10.out
mv image_HNC_uclchem_J10.out uclchem/
mv output_HNC_uclchem_J10.txt uclchem/

#J=2-1
cp -p HNC_J21/camera_wavelength_micron_HNCJ21.inp .
mv camera_wavelength_micron_HNCJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J21.txt
mv image.out image_HNC_uclchem_J21.out
mv image_HNC_uclchem_J21.out uclchem/
mv output_HNC_uclchem_J21.txt uclchem/

#J=3-2
cp -p HNC_J32/camera_wavelength_micron_HNCJ32.inp .
mv camera_wavelength_micron_HNCJ32.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J32.txt
mv image.out image_HNC_uclchem_J32.out
mv image_HNC_uclchem_J32.out uclchem/
mv output_HNC_uclchem_J32.txt uclchem/

: <<'COMMENT'
#J=4-3
cp -p HNC_J43/camera_wavelength_micron_HNCJ43.inp .
mv camera_wavelength_micron_HNCJ43.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J43.txt
mv image.out image_HNC_uclchem_J43.out
mv image_HNC_uclchem_J43.out uclchem/
mv output_HNC_uclchem_J43.txt uclchem/

#J=5-4
cp -p HNC_J54/camera_wavelength_micron_HNCJ54.inp .
mv camera_wavelength_micron_HNCJ54.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J54.txt
mv image.out image_HNC_uclchem_J54.out
mv image_HNC_uclchem_J54.out uclchem/
mv output_HNC_uclchem_J54.txt uclchem/

#J=65
cp -p HNC_J65/camera_wavelength_micron_HNCJ65.inp .
mv camera_wavelength_micron_HNCJ65.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J65.txt
mv image.out image_HNC_uclchem_J65.out
mv image_HNC_uclchem_J65.out uclchem/
mv output_HNC_uclchem_J65.txt uclchem/

#J=76
cp -p HNC_J76/camera_wavelength_micron_HNCJ76.inp .
mv camera_wavelength_micron_HNCJ76.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J76.txt
mv image.out image_HNC_uclchem_J76.out
mv image_HNC_uclchem_J76.out uclchem/
mv output_HNC_uclchem_J76.txt uclchem/

#J=87
cp -p HNC_J87/camera_wavelength_micron_HNCJ87.inp .
mv camera_wavelength_micron_HNCJ87.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J87.txt
mv image.out image_HNC_uclchem_J87.out
mv image_HNC_uclchem_J87.out uclchem/
mv output_HNC_uclchem_J87.txt uclchem/

#J=98
cp -p HNC_J98/camera_wavelength_micron_HNCJ98.inp .
mv camera_wavelength_micron_HNCJ98.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J98.txt
mv image.out image_HNC_uclchem_J98.out
mv image_HNC_uclchem_J98.out uclchem/
mv output_HNC_uclchem_J98.txt uclchem/

#J=109
cp -p HNC_J109/camera_wavelength_micron_HNCJ109.inp .
mv camera_wavelength_micron_HNCJ109.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HNC_uclchem_J109.txt
mv image.out image_HNC_uclchem_J109.out
mv image_HNC_uclchem_J109.out uclchem/
mv output_HNC_uclchem_J109.txt uclchem/
COMMENT

rm radmc3d.out
rm camera_wavelength_micron*
rm numberdens_hnc.inp
rm lines.inp
rm gas_temperature.inp
