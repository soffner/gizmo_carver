#!/bin/bash
#PBS -P iv23
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l ncpus=48
#PBS -l mem=100GB
#PBS -l storage=scratch/jh2+gdata/jh2
#PBS -l wd
#PBS -N carver_dochem
#PBS -j oe
#PBS -m bea
#PBS -M u6645980@alumni.anu.edu.au

python3 do_chemistry_parallel.py > shell.out
