#! /bin/bash

verbose_mode=false
events=0
list_mode=false
listfile="run_list.list"
dosubmit=false
singletest=false
user=`id -u -n`
nruns=0
nseg=0
outdir=$(pwd)"/../data"
run_N=0 #if doing multiple runs at once this gets reset each time


make_condor_jobs() 
{
	#argument 1 = nfiles--derived or given 
	outdir=${outdir}"/run_${run_N}/"
	if [ ! -d ${outdir} ]; then 
		mkdir -p ${outdir}
	fi
	if [ ! -d $(pwd)"/condor_file_dir" ]; then
		mkdir -p $(pwd)"/condor_file_dir"
	fi

	for i in $(seq 0 ${nseg}); do
		j=$(( i + 1 )) 
		condor_file="$(pwd)/condor_file_dir/condor_run-"$run_N"_seg-"$j".job"
		condor_out_file=$(pwd)"/condor_file_dir/condor_run-"$run_N"_seg-"$j".out"
		condor_err_file=$(pwd)"/condor_file_dir/condor_run-"$run_N"_seg-"$j".err"
		condor_log_file="/tmp/Skaydi_condor_run-"$run_N"_seg-"$j".log"
		input_file=`sed "${j}q;d" <data_file_list/dst_calofitting-000${run_N}.list`
		if [ "$vebose_mode" = true ]; then
			echo "Producing condor job file " $condor_file
		fi
		condor_testfile="condor_blank.job"
		IFS=$'\n' read -d '' -r -a blanklines < $condor_testfile
		echo "${blanklines[0]}" > $condor_file 
		echo "${blanklines[1]}"$(pwd)"/RunEMCALChannelTiming.sh" >> $condor_file
		echo "${blanklines[2]}" $input_file $events >> $condor_file
		echo "${blanklines[3]}"$condor_out_file >> $condor_file
		echo "${blanklines[4]}"$condor_err_file >> $condor_file
		echo "${blanklines[5]}"$condor_log_file >> $condor_file
		echo "${blanklines[6]}" $outdir >>$condor_file
		echo "${blanklines[7]}" >> $condor_file
		echo "${blanklines[8]}" >> $condor_file 
		echo "${blanklines[9]}" "   "  $user >> $condor_file 
		echo "${blanklines[10]}" >> $condor_file
		echo "${blanklines[11]}" >> $condor_file
		echo "${blanklines[12]}" >> $condor_file
		echo "${blanklines[13]}" >> $condor_file
	done		
	
}

submit_condor_jobs()
{
	for n in $(seq 0 ${nseg}); do 
		i=$(pwd)"/condor_file_dir/condor_run-"$run_N"_seg-"$n".job"
		condor_submit $i
	done
}

has_argument(){
	[[ ("$1" == *=* && -n ${1#*=}) || ( ! -z "$2" && "$2" != -*) ]]
}

extract_argument() {
	echo "${2:-${1#*=}}"
}

handle_options()
{
	while [ $# -gt 0 ]; do 
		case $1 in 
			-h | --help) 
				echo "Options for EMCAL Chanel Timing Tests"

				echo "$0 [OPTIONS]"
				echo "This script runs the EMCAL Channel Timing Analysis"
				echo ""
				echo " -h, --help 	Display this message"
				echo " -v, --verbose	Enable verbose job creation (Default false) "
				echo " -s, --submit 	Submit condor jobs (default false) "
				echo " -n, --events	Number of events (default 0 runs all) "
				echo " -S, --segments	Number of segments per run (default 0 runs all) " 
				echo " -N, --runs	Number of runs to use (default 0 runs everything from run_list.list)"
				echo " -l, --list	List of runs to use (default run_list.list) "
				echo " -r, --run	Give a specific run number (sets to single run mode) "
				echo " -R, --single	Only run over one run "
				echo " -o, --outdir	Output directory (default ../data/run_[runnumber]) "
				exit 0 
				;;
			-v | --verbose)
				verbose_mode=true
				shift 
				;;
			-s | --submit)
				dosubmit=true
				shift
				;;
			-n | --events)
				events=$(extract_argument $@)
				shift 
				shift
				;;
			-S | --segments)
				nseg=$(extract_argument $@)
				shift
				shift
				;;
			-N | --runs)
				nruns=$(extract_argument $@)
				if [[ $nruns -eq 1 ]]; then
					singlemode=true
				else 
					singlemode=false
				fi

				shift 
				shift
				;;
			-l | --list)
				listfile=$(extract_argument $@)
				listmode=true
				shift
				shift
				;;
			-r | --run)
				nruns=1
				run_N=$(extract_argument $@)
				singlemode=true
				shift
				shift
				;;
			-R | --single)
				nruns=1
				run_N=`head -n 1 ${listfile}`
				singlemode=true
				shift
				shift
				;;
			-o | --outdir)
				outdir=$(extract_argument $@)
				shift 
				shift
				;;
			*)
				echo "Invalid option: $1 "
				exit 1
				::
		esac
	done
}
generate_many_runs()
{
	for i in $(seq 0 ${nruns}); do
		j=$(( i + 1 ))	
		run_N= `sed "${j}q;d" ${listfile}`
		generate_one_run
	done
}
generate_one_run()
{
	hold_nseg=$nseg
	if [[ $nseg -gt `wc -l < data_file_list/dst_calofitting-000${run_N}.list` ]]; then
		nseg=0
	fi
	if [[ $nseg -eq 0 ]]; then
		nseg=`wc -l < data_file_list/dst_calofitting-000${run_N}.list`
	fi
	make_condor_jobs
	if [[ "${dosubmit}" = true ]]; then
		submit_condor_jobs
	fi
	nseg=$hold_nseg #makes sure to not reset between runs 
}
handle_options "$@"
if [[ $nruns -eq 0 ]]; then 
	nruns = `wc -l < ${listfile}`
fi
if [[ "${listmode}" = true ]]; then 
	generate_many_runs
else
	generate_one_run
fi

