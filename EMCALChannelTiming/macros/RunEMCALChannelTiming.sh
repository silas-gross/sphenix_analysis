#! /bin/bash

inputfile=${1:-""}
nevents=${2:-"0"}
MYINSTALL="/sphenix/user/sgross/install_dir"

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana 
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/root-6.32.06/bin/root -x $(pwd)/../../macros/RunEMCALChannelTiming.C\(\"$inputfile\",$nevents\)
