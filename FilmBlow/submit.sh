#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=1000m 
#SBATCH --time=10-02:00:00
#SBATCH --account=rlarson0
#SBATCH --partition=standard

# The application(s) to execute along with its input arguments and options:

srun -n 1 ./FILMBLOW>log.txt
