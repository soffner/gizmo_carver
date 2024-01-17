#!/bin/bash
#PBS -P iv23
#PBS -q normal
#PBS -l walltime=14:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd
#PBS -N carver_HCOp
#PBS -j oe
#PBS -m bea
#PBS -M u6645980@alumni.anu.edu.au

#to write levelpop_co.dat, include writepop keyword: radmc3d image npix 256 loadlambda fluxcons doppcatch inclline linelist nostar writepop sizepc 5.0 phi 0 incl 0 | tee output.txt
sizepc=5.0
phi=0
incl=0
npix=256

#get lines.inp file for HCOp
cp -p lines/lines_hcop.inp .
mv lines_hcop.inp lines.inp

#get temperature data for HCOp
cp -p temperatures/gas_temperature_HCOp.inp .
mv gas_temperature_HCOp.inp gas_temperature.inp

#HCOp Despotic
cp -p despotic/numberdens_HCOp_despotic.inp .
mv numberdens_HCOp_despotic.inp numberdens_hcop.inp

#J=1-0
cp -p HCOp_J10/camera_wavelength_micron_HCOpJ10.inp .
mv camera_wavelength_micron_HCOpJ10.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J10.txt
mv image.out image_HCOp_despotic_J10.out
mv image_HCOp_despotic_J10.out despotic/
mv output_HCOp_despotic_J10.txt despotic/

#J=2-1
cp -p HCOp_J21/camera_wavelength_micron_HCOpJ21.inp .
mv camera_wavelength_micron_HCOpJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J21.txt
mv image.out image_HCOp_despotic_J21.out
mv image_HCOp_despotic_J21.out despotic/
mv output_HCOp_despotic_J21.txt despotic/

#J=3-2
cp -p HCOp_J32/camera_wavelength_micron_HCOpJ32.inp .
mv camera_wavelength_micron_HCOpJ32.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J32.txt
mv image.out image_HCOp_despotic_J32.out
mv image_HCOp_despotic_J32.out despotic/
mv output_HCOp_despotic_J32.txt despotic/

: <<'COMMENT'
#J=4-3
cp -p HCOp_J43/camera_wavelength_micron_HCOpJ43.inp .
mv camera_wavelength_micron_HCOpJ43.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J43.txt
mv image.out image_HCOp_despotic_J43.out
mv image_HCOp_despotic_J43.out despotic/
mv output_HCOp_despotic_J43.txt despotic/

#J=5-4
cp -p HCOp_J54/camera_wavelength_micron_HCOpJ54.inp .
mv camera_wavelength_micron_HCOpJ54.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J54.txt
mv image.out image_HCOp_despotic_J54.out
mv image_HCOp_despotic_J54.out despotic/
mv output_HCOp_despotic_J54.txt despotic/

#J=65
cp -p HCOp_J65/camera_wavelength_micron_HCOpJ65.inp .
mv camera_wavelength_micron_HCOpJ65.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J65.txt
mv image.out image_HCOp_despotic_J65.out
mv image_HCOp_despotic_J65.out despotic/
mv output_HCOp_despotic_J65.txt despotic/

#J=76
cp -p HCOp_J76/camera_wavelength_micron_HCOpJ76.inp .
mv camera_wavelength_micron_HCOpJ76.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J76.txt
mv image.out image_HCOp_despotic_J76.out
mv image_HCOp_despotic_J76.out despotic/
mv output_HCOp_despotic_J76.txt despotic/

#J=87
cp -p HCOp_J87/camera_wavelength_micron_HCOpJ87.inp .
mv camera_wavelength_micron_HCOpJ87.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J87.txt
mv image.out image_HCOp_despotic_J87.out
mv image_HCOp_despotic_J87.out despotic/
mv output_HCOp_despotic_J87.txt despotic/

#J=98
cp -p HCOp_J98/camera_wavelength_micron_HCOpJ98.inp .
mv camera_wavelength_micron_HCOpJ98.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J98.txt
mv image.out image_HCOp_despotic_J98.out
mv image_HCOp_despotic_J98.out despotic/
mv output_HCOp_despotic_J98.txt despotic/

#J=109
cp -p HCOp_J109/camera_wavelength_micron_HCOpJ109.inp .
mv camera_wavelength_micron_HCOpJ109.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_despotic_J109.txt
mv image.out image_HCOp_despotic_J109.out
mv image_HCOp_despotic_J109.out despotic/
mv output_HCOp_despotic_J109.txt despotic/
COMMENT

#HCOp UCLCHEM
cp -p uclchem/numberdens_HCOp_uclchem.inp .
mv numberdens_HCOp_uclchem.inp numberdens_hcop.inp

#J=1-0
cp -p HCOp_J10/camera_wavelength_micron_HCOpJ10.inp .
mv camera_wavelength_micron_HCOpJ10.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl setthreads $PBS_NCPUS | tee output_HCOp_uclchem_J10.txt
mv image.out image_HCOp_uclchem_J10.out
mv image_HCOp_uclchem_J10.out uclchem/
mv output_HCOp_uclchem_J10.txt uclchem/

#J=2-1
cp -p HCOp_J21/camera_wavelength_micron_HCOpJ21.inp .
mv camera_wavelength_micron_HCOpJ21.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J21.txt
mv image.out image_HCOp_uclchem_J21.out
mv image_HCOp_uclchem_J21.out uclchem/
mv output_HCOp_uclchem_J21.txt uclchem/

#J=3-2
cp -p HCOp_J32/camera_wavelength_micron_HCOpJ32.inp .
mv camera_wavelength_micron_HCOpJ32.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J32.txt
mv image.out image_HCOp_uclchem_J32.out
mv image_HCOp_uclchem_J32.out uclchem/
mv output_HCOp_uclchem_J32.txt uclchem/

: <<'COMMENT'
#J=4-3
cp -p HCOp_J43/camera_wavelength_micron_HCOpJ43.inp .
mv camera_wavelength_micron_HCOpJ43.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J43.txt
mv image.out image_HCOp_uclchem_J43.out
mv image_HCOp_uclchem_J43.out uclchem/
mv output_HCOp_uclchem_J43.txt uclchem/

#J=5-4
cp -p HCOp_J54/camera_wavelength_micron_HCOpJ54.inp .
mv camera_wavelength_micron_HCOpJ54.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J54.txt
mv image.out image_HCOp_uclchem_J54.out
mv image_HCOp_uclchem_J54.out uclchem/
mv output_HCOp_uclchem_J54.txt uclchem/

#J=65
cp -p HCOp_J65/camera_wavelength_micron_HCOpJ65.inp .
mv camera_wavelength_micron_HCOpJ65.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J65.txt
mv image.out image_HCOp_uclchem_J65.out
mv image_HCOp_uclchem_J65.out uclchem/
mv output_HCOp_uclchem_J65.txt uclchem/

#J=76
cp -p HCOp_J76/camera_wavelength_micron_HCOpJ76.inp .
mv camera_wavelength_micron_HCOpJ76.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J76.txt
mv image.out image_HCOp_uclchem_J76.out
mv image_HCOp_uclchem_J76.out uclchem/
mv output_HCOp_uclchem_J76.txt uclchem/

#J=87
cp -p HCOp_J87/camera_wavelength_micron_HCOpJ87.inp .
mv camera_wavelength_micron_HCOpJ87.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J87.txt
mv image.out image_HCOp_uclchem_J87.out
mv image_HCOp_uclchem_J87.out uclchem/
mv output_HCOp_uclchem_J87.txt uclchem/

#J=98
cp -p HCOp_J98/camera_wavelength_micron_HCOpJ98.inp .
mv camera_wavelength_micron_HCOpJ98.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J98.txt
mv image.out image_HCOp_uclchem_J98.out
mv image_HCOp_uclchem_J98.out uclchem/
mv output_HCOp_uclchem_J98.txt uclchem/

#J=109
cp -p HCOp_J109/camera_wavelength_micron_HCOpJ109.inp .
mv camera_wavelength_micron_HCOpJ109.inp camera_wavelength_micron.inp
radmc3d image npix $npix loadlambda fluxcons doppcatch inclline linelist nostar sizepc $sizepc phi $phi incl $incl | tee output_HCOp_uclchem_J109.txt
mv image.out image_HCOp_uclchem_J109.out
mv image_HCOp_uclchem_J109.out uclchem/
mv output_HCOp_uclchem_J109.txt uclchem/
COMMENT

rm radmc3d.out
rm camera_wavelength_micron*
rm numberdens_hcop.inp
rm lines.inp
rm gas_temperature.inp
