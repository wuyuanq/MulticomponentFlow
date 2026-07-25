
#!/bin/sh
#@ job_name         = MulticomponentFlow_hpc
#@ output           = $(job_name).$(jobid).out
#@ error            = $(job_name).$(jobid).err
#@ job_type         = parallel
#@ environment      = COPY_ALL
#@ wall_clock_limit = 03:00:00
#@ node             = 1
#@ tasks_per_node   = 4
#@ queue

source /etc/profile.d/modules.sh
module load openmpi

$MPIEXEC -np $LOADL_TOTAL_TASKS ./MulticomponentFlow_hpc
