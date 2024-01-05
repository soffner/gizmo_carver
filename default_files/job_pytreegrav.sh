#!/bin/bash
#PBS -P iv23
#PBS -q express
#PBS -l walltime=24:00:00
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd
#PBS -N carver_pytreegrav_coldens
#PBS -j oe
#PBS -m bea
#PBS -M u6645980@alumni.anu.edu.au

python3 get_coldens_from_pytreegrav.py
