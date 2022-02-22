#!/bin/bash

BASE_OUT_DIR=$1
N_TO_RUN=$2
I_JOB=$3

echo "RUNNING " $BASE_OUT_DIR $N_TO_RUN $I_JOB

# Check if this solves job crashing problem. Typically run 100 jobs, so at most wait 15 minutes before starting.
sleep $(( I_JOB*9 ))s

# Set up SND environment
source /cvmfs/sndlhc.cern.ch/latest/setUp.sh
cd /eos/home-c/cvilela/SND_FEB_18/

eval `alienv load sndsw/latest-feature-muon_reco-release`
cd -

xrdcp /eos/experiment/sndlhc/convertedData/commissioning-h6-NewCalib//sndsw_raw_000284.root .

export EOSSHIP=root://eosuser.cern.ch/
# Run muon reconstruction
python ${SNDSW_ROOT}/shipLHC/run_muonRecoSND.py -f sndsw_raw_000284.root -g /eos/experiment/sndlhc/convertedData/commissioning-h6-NewCalib/geofile_sndlhc_H6.root -p none -n 2000
mv sndsw_raw_000284_muonReco.root sndsw_raw_000284_muonReco_allDetectors.root

python ${SNDSW_ROOT}/shipLHC/run_muonRecoSND.py -f sndsw_raw_000284.root -g /eos/experiment/sndlhc/convertedData/commissioning-h6-NewCalib/geofile_sndlhc_H6.root -p none -n 2000 --no-use_mufi
mv sndsw_raw_000284_muonReco.root sndsw_raw_000284_muonReco_scifiOnly.root

rm -rf sndsw_raw_000284.root

# Copy output
mkdir -p ${BASE_OUT_DIR}/${I_JOB}/
xrdcp ./*.root ${BASE_OUT_DIR}/${I_JOB}/
xrdcp ./*.status ${BASE_OUT_DIR}/${I_JOB}/
rm -rf ./*.root ./*.status
