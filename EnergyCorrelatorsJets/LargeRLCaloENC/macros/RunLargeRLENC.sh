#! /bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n 
MYINSTALL=/sphenix/user/sgross/install_dir
echo $MYINSTALL
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

data_dst=${1:-'none'}
data_fitting_dst=${2:-'none'}
truth_jet_dst=${3:-'none'}
calo_reco_dst=${4:-'none'}
truth_reco_dst=${5:-'none'}
global_reco_dst=${6:-'none'}
runs=${7:-'0'}
minpt=${8:-'1.0'}

/cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/utils/bin/valgrind  --track-origins=yes --log-file=/gpfs/mnt/gpfs02/sphenix/user/sgross/sphenix_analysis/EnergyCorrelatorsJets/LargeRLCaloENC/macros/data_valdrind-out.txt /cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/root-6.26.06.p01/bin/root.exe -x /sphenix/user/sgross/sphenix_analysis/EnergyCorrelatorsJets/LargeRLCaloENC/macros/RunLargeRLENC.C\(\"$data_dst\",\"$data_fitting_dst\",\"$truth_jet_dst\",\"$calo_reco_dst\",\"$truth_reco_dst\",\"$global_reco_dst\",\"$runs\",\"$minpt\"\)

