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
export GALGCONF=/eos/home-c/cvilela/genie_conf/
#eval `alienv load sndsw/latest-96e09ffc70-release`
eval `alienv load sndsw/latest-feature-muon_reco-release`
cd -

# Run detector simulation
python $SNDSW_ROOT/shipLHC/run_simSND.py --Ntuple -n $N_TO_RUN -i $(( I_JOB*N_TO_RUN )) -f /eos/experiment/sndlhc/MonteCarlo/FLUKA/muons_down/muons_VCdown_IR1-LHC.root

export EOSSHIP=root://eosuser.cern.ch/
# Run digitization
python $SNDSW_ROOT/shipLHC/run_digiSND.py -g ./geofile_full.Ntuple-TGeant4.root -f ./sndLHC.Ntuple-TGeant4.root

# Run muon reconstruction
python ${SNDSW_ROOT}/shipLHC/run_muonRecoSND.py -f ./sndLHC.Ntuple-TGeant4_dig.root -p ./ship.params.Ntuple-TGeant4.root -g ./geofile_full.Ntuple-TGeant4.root
mv sndLHC.Ntuple-TGeant4_dig_muonReco.root sndLHC.Ntuple-TGeant4_dig_muonReco_allDetectors.root

python ${SNDSW_ROOT}/shipLHC/run_muonRecoSND.py -f ./sndLHC.Ntuple-TGeant4_dig.root -p ./ship.params.Ntuple-TGeant4.root -g ./geofile_full.Ntuple-TGeant4.root --no-use_mufi
mv sndLHC.Ntuple-TGeant4_dig_muonReco.root sndLHC.Ntuple-TGeant4_dig_muonReco_scifiOnly.root

# Copy output
mkdir -p ${BASE_OUT_DIR}/${I_JOB}/
xrdcp ./*.root ${BASE_OUT_DIR}/${I_JOB}/
xrdcp ./*.status ${BASE_OUT_DIR}/${I_JOB}/
rm -rf ./*.root ./*.status
