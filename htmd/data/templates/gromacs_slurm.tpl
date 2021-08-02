#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --account="{{ account }}"
#SBATCH --partition="{{ partition }}"
#SBATCH --nodes="{{ nodes }}"
#SBATCH --ntasks-per-node="{{ taskspernode }}"
#SBATCH --time="{{ time }}"
#SBATCH --mem="{{ mem }}"

### solver =

#### gromacs run commands here