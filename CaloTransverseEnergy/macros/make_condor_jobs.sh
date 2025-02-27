#! /bin/bash 
#input is run number submit  
RUN=${1:-'21518'}
NFILES=${2:-'0'}
SUBMIT=${3:-'test'}
i=0
LOC=/sphenix/lustre01/sphnxpro/commissioning/DST_ana395_2023p007
FN=$(ls $LOC"/DST_"*$RUN*".root" -l | wc -l)
if [[ $FN -eq 0 ]]; then 
	LOC=/sphenix/lustre01/sphnxpro/commissioning/DST_ana391_2023p006/
fi
for infile in  $LOC/DST_*$RUN*.root; do 
	if [[ $NFILES -ne 0 ]];  then 
		if [[ $i -gt $NFILES ]]; then 
			continue
		fi
	fi
	fname="condor_"$i"_run_"$RUN".job"
	touch $fname
	echo $fname
	
echo "Universe        = vanilla" > $fname
echo "Executable      = /gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/CaloTransverseEnergy/macros/GetET.sh" >> $fname
#	if ((i<10)); then 
#		echo "Arguments       = /sphenix/tg/tg01/jets/ahodges/run23_production/$RUN/DST-000$RUN-000$i.root " >> $fname 
#	else 
	echo "Arguments       = $infile " >> $fname
#	fi
echo "Output          = /gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/CaloTransverseEnergy/running_dir/condor_"$i"_run_"$RUN".out" >> $fname
echo "Error           =/gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/CaloTransverseEnergy/running_dir/condor_"$i"_run_"$RUN".err" >> $fname
echo "Log             =/gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/CaloTransverseEnergy/running_dir/condor_"$i"_run_"$RUN".log" >> $fname
echo "Initialdir      = /gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/CaloTransverseEnergy/src" >> $fname
echo "PeriodicHold    = (NumJobStarts>=1 && JobStatus == 1)" >> $fname
echo "accounting_group = group_phenix.u" >> $fname
echo "accounting_group_user = sgross" >> $fname
echo "request_memory = 128192MB" >> $fname
echo "Priority = 90" >> $fname
echo "job_lease_duration = 3600" >> $fname
echo "Queue 1 " >> $fname


if [[ $SUBMIT == "submit" ]]; then 
	condor_submit $fname
fi
i=$(( i+1 ))
done
 	
