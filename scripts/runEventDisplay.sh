#!/bin/bash

FLAVOUR=${1}
FILE_NUMBER=${2}
EVENT_NUMBER=${3}

BASE_MC_LOCATION=/eos/home-c/cvilela/SND_ANALYSIS/neutrino/

BASE_CODE_LOCATION=/eos/home-c/cvilela/SND_NOV_2/snd_muon_reco/

python ${BASE_CODE_LOCATION}/event_display.py ${BASE_MC_LOCATION}/${FLAVOUR}/${FILE_NUMBER}/sndLHC.Genie-TGeant4_dig.root ${BASE_MC_LOCATION}/${FLAVOUR}/${FILE_NUMBER}/geofile_full.Genie-TGeant4.root $EVENT_NUMBER --reco_file ${BASE_MC_LOCATION}/${FLAVOUR}/${FILE_NUMBER}/sndLHC.Genie-TGeant4_dig_muonReco.root
