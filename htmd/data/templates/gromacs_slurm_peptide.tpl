#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --account="{{ account }}"
#SBATCH --partition="{{ partition }}"
##SBATCH --nodes=1
#SBATCH --ntasks=8
##SBATCH --ntasks-per-node={{ taskspernode }}
#SBATCH --time={{ time }}
#SBATCH --mem={{ mem }}

### will likely need to specify the complete path to each of these inputs

module load icc_18 icc_18-impi_2018 gcc/6.3.1
source /gscratch/pfaendtner/jpfaendt/codes/plumed2/sourceme.sh
source /gscratch/pfaendtner/jpfaendt/codes/gmx2020.5-cpu/bin/GMXRC

##### Energy Minimization #####
gmx_mpi grompp -f {{ em_mdp }} -c {{ system }} -p {{ topol }} -o {{ name }}_em.tpr -maxwarn 3
mpirun -np {{ cores }} gmx_mpi mdrun -ntomp 1 -s {{ name }}_em.tpr -c {{ name }}_em.gro -e {{ name }}_em.edr -g {{ name }}_em.log

##### NPT Equilibration #####
gmx_mpi grompp -f {{ npt_mdp }} -c {{ name }}_em.gro  -p {{ topol }} -o {{ name }}_npt.tpr -maxwarn 3
### when running npt for the first time
mpirun -np {{ cores }} gmx_mpi mdrun -ntomp 1 -s {{ name }}_npt.tpr -nsteps 10 -c {{ name }}_npt.gro -e {{ name }}_npt.edr -x {{ name }}_npt.xtc -g {{ name }}_npt.log -cpi  -cpo  -cpt 1.0 -v >& {{ name }}_log.txt
### once the above is finished
mpirun -np {{ cores }} gmx_mpi mdrun -ntomp 1 -s {{ name }}_npt.tpr -c {{ name }}_npt.gro -e {{ name }}_npt.edr -x {{ name }}_npt.xtc -g {{ name }}_npt.log -cpi  -cpo  -cpt 1.0 -append -v >& {{ name }}_log.txt

##### NVT Equilibration #####
gmx_mpi grompp -f {{ nvt_mdp }} -c {{ name }}_npt.gro -p {{ topol }} -o {{ name }}_nvt.tpr -maxwarn 3
### when running nvt for the first time
mpirun -np {{ cores }} gmx_mpi mdrun -ntomp 1 -s {{ name }}_nvt.tpr -nsteps 10 -c {{ name }}_nvt.gro -e {{ name }}_nvt.edr -x {{ name }}_nvt.xtc -g {{ name }}_nvt.log -cpi restart -cpo restart -cpt 1.0 -v >& {{ name }}_log.txt
### once the obove is finished
mpirun -np {{ cores }} gmx_mpi mdrun -ntomp 1 -s {{ name }}_nvt.tpr -c {{ name }}_nvt.gro -e {{ name }}_nvt.edr -x {{ name }}_nvt.xtc -g {{ name }}_nvt.log -cpi  restart -cpo restart -append -cpt 1.0 -v >& name_log.txt

exit 0
