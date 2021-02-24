#!/bin/bash

cd $SLURM_SUBMIT_DIR

PrjName="deTransAsm"

if [[ -z "${param// }" ]];
then
	Slrm_J="-J $PrjName"
	echo "Setting an slurm job submission name option: $Slrm_J"
else
	Slrm_J=""
	echo "No slurm job submission name set! PrjName=$PrjName"
fi

module load snakemake/5.13.0
module load python/3.7

R=/data/NCBR/projects/SRasmbly/denovo-transcriptome-assembly

modtime1=`stat -c %y $R/Reports/snakemake.log|awk -F "." '{print $1}'|sed 's/ /_/g' -|sed 's/:/_/g'|sed 's/-/_/g' -`
modtime2=`stat -c %y $R/Reports/makeasnake.log|awk -F "." '{print $1}'|sed 's/ /_/g' -|sed 's/:/_/g'|sed 's/-/_/g' -`
mkdir -p $R/Reports

touch $R/Reports/snakemake.log
touch $R/Reports/makeasnake.log

mv $R/Reports/snakemake.log $R/Reports/snakemake.log.$modtime1
mv $R/Reports/makeasnake.log $R/Reports/makeasnake.log.$modtime2

touch $R/Reports/makeasnake.log
touch $R/Reports/snakemake.log

snakemake --unlock -j 1
sbatch $Slrm_J --partition=norm --gres=lscratch:500 --time=10-00:00:00 --mail-type=BEGIN,END,FAIL $R/pipeline_ctrl.sh
