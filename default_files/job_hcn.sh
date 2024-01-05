#!/bin/bash
#PBS -P iv23
#PBS -q express
#PBS -l walltime=18:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd
#PBS -N carver_HCN
#PBS -j oe
#PBS -m bea
#PBS -M u6645980@alumni.anu.edu.au

#to write levelpop_co.dat, include writepop keyword: radmc3d image npix 256 loadlambda fluxcons doppcatch inclline linelist nostar writepop sizepc 5.0 phi 0 incl 0 | tee output.txt
sizepc=5.0
phi=0
incl=0
npix=256

#get lines.inp file for HCN
cp -p lines/lines_hcn.inp .
mv lines_hcn.inp lines.inp

#get temperature data for HCN
cp -p temperatures/gas_temperature_HCN.inp .
mv gas_temperature_HCN.inp gas_temperature.inp

#HCN UCLCHEM
cp -p uclchem/numberdens_HCN_uclchem.inp .
mv numberdens_HCN_uclchem.inp numberdens_hcn.inp

#J=1-0
cp -p HCN_J10/camera_wavelength_micron_HCNJ10.inp .
mv camera_wavelength_micron_HCNJ10.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl setthreads $PBS_NCPUS | tee output_HCN_uclchem_J10.txt
mv image.out image_HCN_uclchem_J10.out
mv image_HCN_uclchem_J10.out uclchem/
mv output_HCN_uclchem_J10.txt uclchem/

#J=2-1
cp -p HCN_J21/camera_wavelength_micron_HCNJ21.inp .
mv camera_wavelength_micron_HCNJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J21.txt
mv image.out image_HCN_uclchem_J21.out
mv image_HCN_uclchem_J21.out uclchem/
mv output_HCN_uclchem_J21.txt uclchem/

#J=3-2
cp -p HCN_J32/camera_wavelength_micron_HCNJ32.inp .
mv camera_wavelength_micron_HCNJ32.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J32.txt
mv image.out image_HCN_uclchem_J32.out
mv image_HCN_uclchem_J32.out uclchem/
mv output_HCN_uclchem_J32.txt uclchem/

: <<'COMMENT'
#J=4-3
cp -p HCN_J43/camera_wavelength_micron_HCNJ43.inp .
mv camera_wavelength_micron_HCNJ43.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J43.txt
mv image.out image_HCN_uclchem_J43.out
mv image_HCN_uclchem_J43.out uclchem/
mv output_HCN_uclchem_J43.txt uclchem/

#J=5-4
cp -p HCN_J54/camera_wavelength_micron_HCNJ54.inp .
mv camera_wavelength_micron_HCNJ54.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J54.txt
mv image.out image_HCN_uclchem_J54.out
mv image_HCN_uclchem_J54.out uclchem/
mv output_HCN_uclchem_J54.txt uclchem/

#J=65
cp -p HCN_J65/camera_wavelength_micron_HCNJ65.inp .
mv camera_wavelength_micron_HCNJ65.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J65.txt
mv image.out image_HCN_uclchem_J65.out
mv image_HCN_uclchem_J65.out uclchem/
mv output_HCN_uclchem_J65.txt uclchem/

#J=76
cp -p HCN_J76/camera_wavelength_micron_HCNJ76.inp .
mv camera_wavelength_micron_HCNJ76.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J76.txt
mv image.out image_HCN_uclchem_J76.out
mv image_HCN_uclchem_J76.out uclchem/
mv output_HCN_uclchem_J76.txt uclchem/

#J=87
cp -p HCN_J87/camera_wavelength_micron_HCNJ87.inp .
mv camera_wavelength_micron_HCNJ87.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J87.txt
mv image.out image_HCN_uclchem_J87.out
mv image_HCN_uclchem_J87.out uclchem/
mv output_HCN_uclchem_J87.txt uclchem/

#J=98
cp -p HCN_J98/camera_wavelength_micron_HCNJ98.inp .
mv camera_wavelength_micron_HCNJ98.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J98.txt
mv image.out image_HCN_uclchem_J98.out
mv image_HCN_uclchem_J98.out uclchem/
mv output_HCN_uclchem_J98.txt uclchem/

#J=109
cp -p HCN_J109/camera_wavelength_micron_HCNJ109.inp .
mv camera_wavelength_micron_HCNJ109.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCN_uclchem_J109.txt
mv image.out image_HCN_uclchem_J109.out
mv image_HCN_uclchem_J109.out uclchem/
mv output_HCN_uclchem_J109.txt uclchem/
COMMENT

rm radmc3d.out
rm camera_wavelength_micron*
rm numberdens_hcn.inp
rm lines.inp
rm gas_temperature.inp
