#! /bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${LOGNAME}/macros/detectors/sPHENIX/

executable="RunEtaShiftStudy.sh"
ex_tag="EtaShift"
nfiles=0
optione=""

make_condor_jobs()
{
	if [[ $nfiles -eq 0 ]]; then 
		nfiles=`wc -l < ${triggertype}_data/jet_density_${filedensity}.list`
	fi	
	for i in $(seq 0 ${nfiles}); do 
		j=$(( i+1 ))
		if [ $i -eq $nfiles ]; then 
			break
		fi
		condor_file="$(pwd)/condor_file_dir/condor_"$triggertype"_seg_"$i".job"
		condor_out_file=$(pwd)"/condor_file_dir/condor_"$triggertype"_seg_"$i".out"
		condor_err_file=$(pwd)"/condor_file_dir/condor_"$triggertype"_seg_"$i".err"
		condor_log_file=$(pwd)"/condor_file_dir/condor_"$triggertype"_seg_"$i".log"
		global=`sed "${j}q;d" ${triggertype}_data/global_density_${filedensity}.list`
		truth=`sed "${j}q;d" ${triggertype}_data/truth_density_${filedensity}.list`
		jet=`sed "${j}q;d" ${triggertype}_data/jet_density_${filedensity}.list`
		calo=`sed "${j}q;d" ${triggertype}_data/calo_density_${filedensity}.list`
		
		if [ "$vebose_mode" = true ]; then
			echo "Producing condor job file " $condor_file
		fi
		IFS=$'\n' read -d '' -r -a blanklines < $condor_testfile
		echo "${blanklines[0]}" > $condor_file 
		echo "${blanklines[1]}"$(pwd)"/${executable}" >> $condor_file
		echo "${blanklines[2]}"${options} >> $condor_file
		echo "${blanklines[3]}"$condor_out_file >> $condor_file
		echo "${blanklines[4]}"$condor_err_file >> $condor_file
		echo "${blanklines[5]}"$condor_log_file >> $condor_file
		echo "${blanklines[6]} $outDir" >>$condor_file
		echo "${blanklines[7]}" >> $condor_file
		echo "${blanklines[8]}" >> $condor_file 
		echo "${blanklines[9]}" "   "  $USER >> $condor_file 
		echo "${blanklines[10]}" >> $condor_file
		echo "${blanklines[11]}" >> $condor_file
		echo "${blanklines[12]}" >> $condor_file
		echo "${blanklines[13]}" >> $condor_file
	done

}

choose_executable() 
{
	#choose the executable 
	if [[ "${ex_tag}" = *"EtaShift"* ]]; then 
		executable="RunEtaShiftStudy.sh"
	fi
}
select_options()
{
	if [[ "${ex_tag}" = *"EtaShift" ]]; then 
		options=""
	fi
}
