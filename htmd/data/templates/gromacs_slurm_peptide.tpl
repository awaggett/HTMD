#!/bin/bash
#SBATCH --job-name="{{ name }}"
#SBATCH --account="{{ account }}"
#SBATCH --partition="{{ partition }}"
#SBATCH --nodes={{ nodes }}
#SBATCH --ntasks-per-node={{ taskspernode }}
#SBATCH --time={{ time }}
#SBATCH --mem={{ mem }}

### will likely need to specify the complete path to each of these inputs

module load icc_18 icc_18-impi_2018 gcc/6.3.1
source /gscratch/pfaendtner/jpfaendt/codes/plumed2/sourceme.sh
source /gscratch/pfaendtner/jpfaendt/codes/gmx2020.5-cpu/bin/GMXRC

name={{ name }}

##### Energy Minimization #####
gmx_mpi grompp -f {{ em_mdp }} -c {{ system }} -p {{ topol }} -o $name_em.tpr -maxwarn 3
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s $name_em.tpr -c $name_em.gro -e $name_em.edr -g $name_em.log }}

##### NPT Equilibration #####
gmx_mpi grompp -f {{ npt_mdp }} -c $name_em.gro  -p {{ topol }} -n {{ index }} -o $name_npt.tpr -maxwarn 3
### when running npt for the first time
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s $name_npt.tpr -nsteps 10 -c $name_npt.gro -e $name_npt.edr -x $name_npt.xtc -g $name_npt.log -cpi  -cpo  -cpt 1.0 -plumed {{ plumed }}  -v >& $name_log.txt
### once the above is finished
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s $name_npt.tpr -c $name_npt.gro }} -e $name_npt.edr -x $name_npt.xtc -g $name_npt.log -cpi  -cpo  -cpt 1.0 -append -plumed {{ plumed }}  -v >& $name_log.txt

##### NVT Equilibration #####
gmx_mpi grompp -f {{ nvt_mdp }} -c $name_npt.gro -p {{ topol }} -n {{ index }} -o $name_nvt.tpr -maxwarn 3
### when running nvt for the first time
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s $name_nvt.tpr -nsteps 10 -c $name_nvt.gro -e $name_nvt.edr -x $name_nvt.xtc -g $name_nvt.log -cpi restart -cpo restart -cpt 1.0 -plumed {{ plumed }} -v >& $name_log.txt
### once the obove is finished
mpirun -np 8 gmx_mpi mdrun -ntomp 1 -s $name_nvt.tpr -c $name_nvt.gro -e $name_nvt.edr -x $name_nvt.xtc -g $name_nvt.log -cpi  restart -cpo restart -append -cpt 1.0 -plumed {{ plumed }} -v >& name_log.txt

exit 0
