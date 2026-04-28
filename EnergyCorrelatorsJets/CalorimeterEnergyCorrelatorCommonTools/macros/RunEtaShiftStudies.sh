#! /bin/bash

thisDir=${1:-""}
caloFile=${2:-""}
globalFile=${3:-""}
jetFile=${4:-""}
truthFile=${5:-""}
outDir=${6:-""}
installdir=${7:-$MYINSTALL}

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.542
source /opt/sphenix/core/bin/setup_local.sh $installdir

n_files=`wc -l < ${caloFile}`
/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/root-6.32.06/bin/root -b -q ${thisDir}/RunEtaShiftStudies.C\(\"$caloFile\",\"$globalFile\",\"$jetFile\",\"$truthFile\",\"$outDir\",\"$n_files\"\)

