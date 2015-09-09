#PBS -S /bin/bash
#PBS -o Fidelity_P93_24.txt
#PBS -e Fidelity_P93_24_err.txt
#PBS -m ae
#PBS -M mwalke49@jhmi.edu
#PBS -q dque
#PBS -l nodes=64:ppn=8
#PBS -l walltime=72:00:00
#PBS -N urgent_job_2

# load the modules
. /etc/profile.modules
module load openmpi/openmpi-1.3.2-intel intel/intel-icc-ifort-Compiler-11.083

MPI_DIR=/apps/openmpi-1.3.2/bin

WORK_DIR=/home/mwalker/models/CRU3D_CPU
EXECUTABLE=./cru3d
PARAM_FILE=param/Fidelity_P93_24.xml
MESH_BASE=meshes/RyR_P93_24

cd $WORK_DIR
make clean;
make;

echo pbs nodefile:
cat $PBS_NODEFILE
NPROCS=`wc -l < $PBS_NODEFILE`
cat $PBS_NODEFILE > host.list

date
time $MPI_DIR/mpirun --mca btl openib,self -np $NPROCS -hostfile host.list $EXECUTABLE -param $PARAM_FILE -mesh $MESH_BASE
date
