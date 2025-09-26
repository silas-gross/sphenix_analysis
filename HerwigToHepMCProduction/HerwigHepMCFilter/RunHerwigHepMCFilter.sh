#! /bin/bash
export MYINSTALL=/sphenix/user/sgross/install_dir
source /opt/sphenix/core/bin/sphenix_setup.sh -n ana
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
input_file=${1:-''}
jet_trigger=${2:-'10.'}
n_events=${3:-'1000'}
root.exe -x -b  $(pwd)/../HerwigHepMCFilter/RunHerwigHepMCFilter.C\(\"$input_file\",\"$jet_trigger\",$n_events\) 
