#! /bin/bash

submit=${1:-'test'}
sample=${2:-'20'}
nFile=${3:-'0'}
cluster=${4:-'false'}
nevts=${5:-'0'}
minpt=${6:-'1.0'}
i=0
if [[ $nFile -eq 0 ]]; then 
	nFile=`wc -l < Pythia${sample}GeVJets/dst_truth_jet.list`
fi 
for i in $(seq 0 ${nFile}); do 
	j=$(( i+1 ))
	outdir=/sphenix/user/sgross/sphenix_analysis/EnergyCorrelatorsJets/LargeRLCaloENC/root_output_Pythia${sample}GeV
        if [ ! -d  $outdir ]; then
                mkdir $outdir
        fi
	fname="condor_files/condor_segment_"$i"_${sample}GeV_cluster_"${cluster}".job"
	data="none" 	 
	truthf="none"
	truthfr=`sed "${j}q;d" Pythia${sample}GeVJets/g4hits.list`
	truthj=`sed "${j}q;d" Pythia${sample}GeVJets/dst_truth_jet.list`
	caloclusterf=`sed "${j}q;d" Pythia${sample}GeVJets/dst_calo_cluster.list`
	globalf=`sed "${j}q;d" Pythia${sample}GeVJets/dst_global.list`
	
	echo "Universe 	        = vanilla " > $fname
	echo "Executable 	= /gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/EnergyCorrelatorsJets/LargeRLCaloENC/macros/RunLargeRLENC.sh " >>$fname
	echo "Arguments         = ${data} none none none ${truthf} ${truthj} ${caloclusterf} ${truthfr} ${globalf} ${cluster} ${nevts} ${minpt}" >> $fname 
	echo "Output  	        = /gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/EnergyCorrelatorsJets/LargeRLCaloENC/macros/condor_files/condor_${i}_${sample}_cluster_${cluster}.out " >> $fname
	echo "Error 		= /gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/EnergyCorrelatorsJets/LargeRLCaloENC/macros/condor_files/condor_${i}_${sample}_cluster_${cluster}.err " >> $fname
	echo "Log  		= /gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/EnergyCorrelatorsJets/LargeRLCaloENC/macros/condor_files/condor_${i}_${sample}_cluster_${cluster}.log" >> $fname
	echo "Initialdir  	= /gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/EnergyCorrelatorsJets/LargeRLCaloENC/root_output_Pythia${sample}GeV" >> $fname
	echo "PeriodicHold 	= (NumJobStarts>=1 && JobStatus == 1)" >>$fname
	echo "accounting_group = group_sphenix.user " >> $fname
	echo "accounting_group_user = sgross " >> $fname
	echo "request_memory = 8 GB " >> $fname
	echo "Priority = 90 ">> $fname
	echo "job_lease_duration = 3600" >> $fname
	echo "Queue 1" >> $fname 
if [[ $submit == "submit" ]]; then 
	condor_submit $fname
fi 		
done 
