#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --account="{{ account }}"
#SBATCH --partition="{{ partition }}"
#SBATCH --nodes="{{ nodes }}"
#SBATCH --ntasks-per-node="{{ taskspernode }}"
#SBATCH --time="{{ time }}"
#SBATCH --mem="{{ mem }}"

### will likely need to specify the complete path to each of these inputs

module load icc_18 icc_18-impi_2018 gcc/6.3.1
source /gscratch/pfaendtner/jpfaendt/codes/plumed2/sourceme.sh
source /gscratch/pfaendtner/jpfaendt/codes/gmx2020.5-cpu/bin/GMXRC

gmx_mpi grompp -f {{ em.mdp }} -c {{ system.gro }} -p {{ topol.top }} -n {{ system.ndx }} -o {{ em.tpr }} -maxwarn 3

mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s {{ em.tpr }} -c {{ em.gro }} -e {{ em.edr }} -g {{ em.log }}

gmx_mpi grompp -f {{ npt.mdp }} -c {{ em.gro }}  -p {{ topol.top }} -n {{ system.ndx }} -o {{ npt.tpr }} -maxwarn 3

### when running npt for the first time
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s {{ npt.tpr }} -nsteps 10 -c {{ npt.gro }} -e {{ npt.edr }} -x {{ npt.xtc }} -g {{ npt.log }} -cpi  -cpo  -cpt 1.0 -plumed {{ plumed.dat }}  -v >& {{ log.txt }}

### once the above is finished
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s {{ npt.tpr }} -c {{ npt.gro }} -e {{ npt.edr }} -x {{ npt.xtc }} -g {{ npt.log }} -cpi  -cpo  -cpt 1.0 -append -plumed {{ plumed.dat }}  -v >& {{ log.txt }}

gmx_mpi grompp -f {{ nvt.mdp }} -c {{ npt.gro }} -p {{ topol.top }} -n {{ system.ndx }} -o {{ nvt.tpr }} -maxwarn 3

### when running nvt for the first time
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s {{ nvt.tpr }} -nsteps 10 -c {{ nvt.gro }} -e {{ nvt.edr }} -x {{ nvt.xtc }} -g {{ nvt.log }} -cpi restart -cpo restart -cpt 1.0 -plumed {{ plumed.dat }} -v >& {{ log.txt }}

### once the obove is finished
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s {{ nvt.tpr }} -c {{ nvt.gro }} -e {{ nvt.edr }} -x {{ nvt.xtc }} -g {{ nvt.log }} -cpi  restart -cpo restart -append -cpt 1.0 -plumed {{ plumed.dat }} -v >& {{ log.txt }}

exit 0
