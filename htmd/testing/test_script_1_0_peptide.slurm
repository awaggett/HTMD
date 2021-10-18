#!/bin/bash
#SBATCH --job-name="test_script_1_0_peptide"
#SBATCH --account="pfaendtner"
#SBATCH --partition="pfaendtner"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=20GB

### will likely need to specify the complete path to each of these inputs

module load icc_18 icc_18-impi_2018 gcc/6.3.1
source /gscratch/pfaendtner/jpfaendt/codes/plumed2/sourceme.sh
source /gscratch/pfaendtner/jpfaendt/codes/gmx2020.5-cpu/bin/GMXRC

##### Energy Minimization #####
gmx_mpi grompp -f /gscratch/pfaendtner/awaggett/htmd_test/HTMD/htmd/data/input_files/em_pep.mdp -c test_script_1_0_peptide_ion.gro -p test_script_1_0_peptide_topol.top -o test_script_1_0_peptide_em.tpr -maxwarn 3
mpirun -np 3 gmx_mpi mdrun -ntomp 1 -s test_script_1_0_peptide_em.tpr -c test_script_1_0_peptide_em.gro -e test_script_1_0_peptide_em.edr -g test_script_1_0_peptide_em.log

##### NPT Equilibration #####
gmx_mpi grompp -f /gscratch/pfaendtner/awaggett/htmd_test/HTMD/htmd/data/input_files/npt_pep.mdp -c test_script_1_0_peptide_em.gro  -p test_script_1_0_peptide_topol.top -o test_script_1_0_peptide_npt.tpr -maxwarn 3
### when running npt for the first time
mpirun -np 3 gmx_mpi mdrun -ntomp 1 -s test_script_1_0_peptide_npt.tpr -nsteps 10 -c test_script_1_0_peptide_npt.gro -e test_script_1_0_peptide_npt.edr -x test_script_1_0_peptide_npt.xtc -g test_script_1_0_peptide_npt.log -cpi  -cpo  -cpt 1.0 -plumed None  -v >& test_script_1_0_peptide_log.txt
### once the above is finished
mpirun -np 3 gmx_mpi mdrun -ntomp 1 -s test_script_1_0_peptide_npt.tpr -c test_script_1_0_peptide_npt.gro -e test_script_1_0_peptide_npt.edr -x test_script_1_0_peptide_npt.xtc -g test_script_1_0_peptide_npt.log -cpi  -cpo  -cpt 1.0 -append -plumed None  -v >& test_script_1_0_peptide_log.txt

##### NVT Equilibration #####
gmx_mpi grompp -f /gscratch/pfaendtner/awaggett/htmd_test/HTMD/htmd/data/input_files/nvt_pep.mdp -c test_script_1_0_peptide_npt.gro -p test_script_1_0_peptide_topol.top -o test_script_1_0_peptide_nvt.tpr -maxwarn 3
### when running nvt for the first time
mpirun -np 3 gmx_mpi mdrun -ntomp 1 -s test_script_1_0_peptide_nvt.tpr -nsteps 10 -c test_script_1_0_peptide_nvt.gro -e test_script_1_0_peptide_nvt.edr -x test_script_1_0_peptide_nvt.xtc -g test_script_1_0_peptide_nvt.log -cpi restart -cpo restart -cpt 1.0 -plumed None -v >& test_script_1_0_peptide_log.txt
### once the obove is finished
mpirun -np 3 gmx_mpi mdrun -ntomp 1 -s test_script_1_0_peptide_nvt.tpr -c test_script_1_0_peptide_nvt.gro -e test_script_1_0_peptide_nvt.edr -x test_script_1_0_peptide_nvt.xtc -g test_script_1_0_peptide_nvt.log -cpi  restart -cpo restart -append -cpt 1.0 -plumed None -v >& name_log.txt

exit 0