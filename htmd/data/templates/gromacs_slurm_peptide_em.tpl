#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --account="{{ account }}"
#SBATCH --partition="{{ partition }}"
#SBATCH --nodes="{{ nodes }}"
#SBATCH --ntasks-per-node="{{ taskspernode }}"
#SBATCH --time="{{ time }}"
#SBATCH --mem="{{ mem }}"

module load icc_18 icc_18-impi_2018 gcc/6.3.1
source /gscratch/pfaendtner/jpfaendt/codes/plumed2/sourceme.sh
source /gscratch/pfaendtner/jpfaendt/codes/gmx2020.5-cpu/bin/GMXRC

##### Energy Minimization #####
##### gmx_mpi grompp -f {{ em.mdp }} -c {{ system.gro }} -p {{ topol.top }} -n {{ system.ndx }} -o {{ em.tpr }} -maxwarn 3
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s {{ em.tpr }} -c {{ em.gro }} -e {{ em.edr }} -g {{ em.log }}

exit 0
