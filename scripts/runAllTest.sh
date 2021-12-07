rm -rf test
mkdir test
export GALGCONF=/eos/home-c/cvilela/SND_ALT/my_genie_conf/
# GENIE
python $SNDSW_ROOT/macro/makeSNDGenieEvents.py sim --nupdg 14 -p CCDIS -n 100 -c /eos/user/c/cvilela/SND_ALT/genie_splines/ -o test/
# DETECTOR SIMULATION
python $SNDSW_ROOT/shipLHC/run_simSND.py --Genie 1 -n 100 -i 0 -f test/numu_CCDIS_FairShip.root -o test/
# DIGITISATION
cd test
python $SNDSW_ROOT/shipLHC/run_digiSND.py -g geofile_full.Genie-TGeant4.root -f sndLHC.Genie-TGeant4.root
cd -
