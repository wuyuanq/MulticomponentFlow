
# @ shell=/bin/csh
# @ job_name = MulticomponentFlow_hpc
# @ account_no = k234
# @ error = error.$(jobid)
# @ output = output.$(jobid)
# @ environment = COPY_ALL;
# @ wall_clock_limit = 1:00:00
# @ notification = always
# @ job_type = bluegene
# @ bg_size = 64
# @ queue

/bgsys/drivers/ppcfloor/bin/mpirun -mode VN -np 4 -exe ./MulticomponentFlow_hpc