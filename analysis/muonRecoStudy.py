import ROOT
import numpy as np
import glob

from collections import OrderedDict

def fiducialVolume(x, y, z) :
    # Target station dimensions. Consider making this more precise.
    if x < -47 :
        return False
    if x > -8 :
        return False
    if y < 15.5 :
        return False
    if y > 54.5 :
        return False
    if z < -25. :
        return False
    if z > 25. :
        return False
    return True

def trueMuons(event) :
    # Define true muon as a muon that starts in the target region and produces at least one hit in the downstream muon filter.
    muons = []

    for hit in event.Digi_MuFilterHits :
        if hit.GetSystem() != 3 :
            continue
        else :
            print("DS MuFilter Hit!")
#            print(event.Digi_MuFilterHits2MCPoints[0].wList(hit.GetDetectorID()))
            for mc_point_i, _ in event.Digi_MuFilterHits2MCPoints[0].wList(hit.GetDetectorID()) :
#                print(mc_point_i)
#                print(event.MuFilterPoint[mc_point_i].PdgCode())
                # Check there's a muon contributing to this hit
                if abs(event.MuFilterPoint[mc_point_i].PdgCode() != 13) :
                    continue
                
                # Check if muon originates in target station
                trackID = event.MuFilterPoint[mc_point_i].GetTrackID()
                track = event.MCTrack[trackID]
                if not fiducialVolume(track.GetStartX(),
                                      track.GetStartY(),
                                      track.GetStartZ()) :
                    continue

                print("Muon ID", trackID)
                if trackID not in muons :
                    muons.append(trackID)

    return muons

# Check if true event vertex is in target station. Also return event weight.
def getNuWeightFV(event) :
    for track in event.MCTrack :
        if track.GetMotherId() == -1 :
            break
    weight = track.GetWeight()
    
    return (fiducialVolume(track.GetStartX(), track.GetStartY(), track.GetStartZ()), weight)

signal_definition_true_cuts = OrderedDict()
event_selection_reco_cuts = OrderedDict()
def cutTracker(dictionary, cut_name, weight = 1) :
    if cut_name not in dictionary :
        dictionary[cut_name] = weight
    else :
        dictionary[cut_name] += weight

from rootpyPickler import Unpickler

sample_name = "NC"
geo_file = "/eos/user/c/cvilela/SND_ANALYSIS/neutrino/"+sample_name+"/*/geofile_full.Genie-TGeant4.root"
input_filename_expr = ["/eos/user/c/cvilela/SND_ANALYSIS/neutrino/"+sample_name+"/*/sndLHC.Genie-TGeant4"]
# Sort out geometry
fgeo = ROOT.TFile.Open(glob.glob(geo_file)[0])
upkl    = Unpickler(fgeo)
snd_geo = upkl.load('ShipGeo')
run = ROOT.FairRunSim()

import shipLHC_conf as sndDet_conf
modules = sndDet_conf.configure(run,snd_geo)
sGeo = fgeo.FAIRGeom
modules['Scifi'].SiPMmapping()

# Sort out input files
input_digi_filenames = []
input_reco_filenames = []
for expr in input_filename_expr :
    for fname in glob.glob(expr+"_dig.root") :
        input_digi_filenames.append(fname)
    for fname in glob.glob(expr+"_dig_muonReco.root") :
        input_reco_filenames.append(fname)

tree_digi = ROOT.TChain("cbmsim")
tree_reco = ROOT.TChain("cbmsim")

for fname in input_digi_filenames :
    print("Adding {0}".format(fname))
    tree_digi.Add(fname)

for fname in input_reco_filenames :
    print("Adding {0}".format(fname))
    tree_reco.Add(fname)

tree_reco.AddFriend(tree_digi)

# Variables for plots
w = []
n_true_mu = []
n_reco_mu = []

# Main event loop
for i_event, event in enumerate(tree_reco) :
    print("Looping through event {0}".format(i_event))
    
    inFV, weight = getNuWeightFV(event)
    cutTracker(signal_definition_true_cuts, "Generated", weight)

    if not inFV :
        print("not in FV", weight)
        continue
    else :
        print("Event is in FV. Weight: {0}".format(weight))
        
    cutTracker(signal_definition_true_cuts, "Vertex in target station", weight)

    true_muons = trueMuons(event)
    print(true_muons)

    # Reco muons!
    n_reco_muons = len(event.Reco_MuonTracks)
    print("Number of reconstructed muons {0}".format(n_reco_muons))

    n_true_mu.append(len(true_muons))
    n_reco_mu.append(n_reco_muons)
    w.append(weight)

    if n_true_mu[-1] == 0 and n_reco_mu[-1] > 0 :
        print("File {0} event {1}".format(event.GetFile().GetPath(), event.GetEntry()))

    # At least one true muon in the event
    if len(true_muons) == 0 :
        continue

    cutTracker(signal_definition_true_cuts, "At least one true muon", weight)
    cutTracker(event_selection_reco_cuts, "At least one true muon", weight)
    
    if n_reco_muons == 0 :
        continue
        
    cutTracker(event_selection_reco_cuts, "At least one reco muon", weight)

# Print some numbers
n_true_mu = np.array(n_true_mu)
n_reco_mu = np.array(n_reco_mu)
w = np.array(w)

for n_mu in range(0, 5) :
    print("True {0} muon events {1} correctly reconstructed {2}".format(n_mu, np.sum(w[n_true_mu == n_mu]), np.sum(w[np.logical_and(n_true_mu == n_mu, n_reco_mu == n_mu)])))
print("Total events {0}".format(np.sum(w)))

# Make some plots
import matplotlib.pyplot as plt

plt.figure()
plt.plot(range(len(signal_definition_true_cuts)), [ w/next(iter(signal_definition_true_cuts.values())) for w in signal_definition_true_cuts.values()])
plt.xticks(range(len(signal_definition_true_cuts)), signal_definition_true_cuts.keys())
plt.tight_layout()

plt.figure()
plt.plot(range(len(event_selection_reco_cuts)), [w/next(iter(event_selection_reco_cuts.values())) for w in event_selection_reco_cuts.values()] )
plt.xticks(range(len(event_selection_reco_cuts)), event_selection_reco_cuts.keys())
plt.tight_layout()

plt.figure()
plt.hist2d(n_true_mu, n_reco_mu, weights = w, range = ((0, 5), (0, 5)), bins = (5, 5))
plt.xlabel("Number of true muons")
plt.ylabel("Number of reconstructed muons")
plt.tight_layout()
plt.savefig("muon_multiplicity_"+sample_name+".png")

plt.show()

