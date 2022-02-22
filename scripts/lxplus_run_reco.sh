#!/bin/bash

BASE_OUT_DIR=$1
N_TO_RUN=$2
I_JOB=$3

NUPDG=$4
NUINT=$5

if [ $NUPDG -eq 12 ]; then NUNAME=nue; elif [ $NUPDG -eq 14 ]; then NUNAME=numu; elif [ $NUPDG -eq -12 ]; then NUNAME=anue; elif [ $NUPDG -eq -14 ]; then NUNAME=anumu;  fi

echo "RUNNING " $BASE_OUT_DIR $N_TO_RUN $I_JOB

# Check if this solves job crashing problem. Typically run 100 jobs, so at most wait 15 minutes before starting.
sleep $(( I_JOB*9 ))s

# Set up SND environment
source /cvmfs/sndlhc.cern.ch/latest/setUp.sh
cd /eos/home-c/cvilela/SND_FEB_7/
export GALGCONF=/eos/home-c/cvilela/genie_conf/
#eval `alienv load sndsw/latest-96e09ffc70-release`
eval `alienv load sndsw/latest-feature-muon_reco-release`
cd -

export EOSSHIP=root://eosuser.cern.ch/

# Run genie
#python $SNDSW_ROOT/macro/makeSNDGenieEvents.py sim --nupdg ${NUPDG} -p ${NUINT} -n $N_TO_RUN -c /eos/home-c/cvilela/genie_splines/ -o ./

# Run detector simulation
#python $SNDSW_ROOT/shipLHC/run_simSND.py --Genie 1 -n ${N_TO_RUN}  -f ./${NUNAME}_${NUINT}_FairShip.root -o ./

# Run digitization
#python $SNDSW_ROOT/shipLHC/run_digiSND.py -g ./geofile_full.Genie-TGeant4.root -f ./sndLHC.Genie-TGeant4.root

xrdcp ${BASE_OUT_DIR}/${I_JOB}/sndLHC.Genie-TGeant4_dig.root .
xrdcp ${BASE_OUT_DIR}/${I_JOB}/ship.params.Genie-TGeant4.root .
xrdcp ${BASE_OUT_DIR}/${I_JOB}/geofile_full.Genie-TGeant4.root .

# Run muon reconstruction
python ${SNDSW_ROOT}/shipLHC/run_muonRecoSND.py -f ./sndLHC.Genie-TGeant4_dig.root -p ./ship.params.Genie-TGeant4.root -g ./geofile_full.Genie-TGeant4.root

# Copy output
rm -rf ${BASE_OUT_DIR}/${I_JOB}/sndLHC.Genie-TGeant4_dig_muonReco.root
xrdcp ./sndLHC.Genie-TGeant4_dig_muonReco.root ${BASE_OUT_DIR}/${I_JOB}/

# Clean up
rm -rf ./*.root ./*.status
