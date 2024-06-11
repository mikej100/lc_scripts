#!/bin/bash
#SBATCH --partition=short
###SBATCH --job-name=SRA
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=00-00:01:00
#SBATCH --output=slurm/%j_%x.out 
#SBATCH --error=slurm/%j_%x.out
##SBATCH --mail-type=fail
#SBATCH --export=TEST_VAR="test value"
#
# Temporary script for feature development.
#
echo "Test environment variable value: ${TEST_VAR}"

now() {
    date +"%Y-%m-%dT%T"
}

echo "$(now) Starting $(basename ${BASH_SOURCE})"

# Show git info for scripts folder
${SCRIPTS}/scripts_info.sh || true

echo "$(now) Ending $(basename ${BASH_SOURCE})"