import ROOT
import numpy as np
import glob

from collections import OrderedDict

# To store detector positions
a = ROOT.TVector3()
b = ROOT.TVector3()

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

sample_name = "numuDefault"
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
plot_vars = {}
plot_vars["w"] = []
plot_vars["n_true_mu"] = []
plot_vars["n_reco_mu_raw"] = []
plot_vars["n_reco_mu_pruned"] = []
plot_vars["n_unassigned_hits"] = []
plot_vars["n_ds_hits"] = []
plot_vars["n_last_ds_hits"] = []
plot_vars["transverse_distance_reco_mu_shower_center"] = []

# Misreco
misid_filename = []
misid_evtno = []
misid_ntrue = []
misid_nreco = []

# Main event loop
for i_event, event in enumerate(tree_reco) :
    if i_event > 20000 :
        break
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
    n_raw_reco_mu = len(event.Reco_MuonTracks)
    print("Number of reconstructed muons {0}".format(n_raw_reco_mu))

    plot_vars["n_true_mu"].append(len(true_muons))
    plot_vars["n_reco_mu_raw"].append(n_raw_reco_mu)
    plot_vars["w"].append(weight)

    # Estimate average position using scifi hits
    scifi_mean_x = 0.
    scifi_mean_y = 0.
    scifi_mean_z = 0.
    
    n_scifi_hor = 0
    n_scifi_ver = 0

    for scifi_hit in event.Digi_ScifiHits :
        modules["Scifi"].GetSiPMPosition(scifi_hit.GetDetectorID(), a, b)
        if scifi_hit.isVertical() :
            n_scifi_ver += 1
            scifi_mean_z += a.Z()
            scifi_mean_x += a.X()
        else :
            n_scifi_hor += 1
            scifi_mean_z += a.Z()
            scifi_mean_y += a.Y()
        
    if n_scifi_ver :
        scifi_mean_x /= n_scifi_ver
    if n_scifi_hor :
        scifi_mean_y /= n_scifi_hor
    if (n_scifi_ver+n_scifi_hor) :
        scifi_mean_z /= (n_scifi_ver+n_scifi_hor)

    # Calculate distance from muon track to shower mean
    track_scifi_d = []
    n_prune_reco_mu = 0
    for track in event.Reco_MuonTracks :
        track_start = track.getStart()
        track_stop = track.getStop()
        track_slope_x = (track_stop.X() - track_start.X())/(track_stop.Z() - track_start.Z())
        track_intercept_x = track_start.X() - track_slope_x*track_start.Z()
        track_slope_y = (track_stop.Y() - track_start.Y())/(track_stop.Z() - track_start.Z())
        track_intercept_y = track_start.Y() - track_slope_y*track_start.Z()
        
        dx = track_slope_x*scifi_mean_z + track_intercept_x - scifi_mean_x
        dy = track_slope_x*scifi_mean_z + track_intercept_y - scifi_mean_y
        
        d = (dx**2 + dy**2)**0.5

        track_scifi_d.append(d)
        
        if d < 10 :
            n_prune_reco_mu += 1
    plot_vars["n_reco_mu_pruned"].append(n_prune_reco_mu)

    if len(track_scifi_d) :
        plot_vars["transverse_distance_reco_mu_shower_center"].append(max(track_scifi_d))
    else :
        plot_vars["transverse_distance_reco_mu_shower_center"].append(0)

    # Count unassigned downstream hits
    n_ds_hits_unassigned = 0
    n_ds_hits = 0
    n_ds_hits_last = 0
    hit_list = []
    for track in event.Reco_MuonTracks :
        for hit in track.getHits() :
            hit_list.append(hit)

    for mufilter_hit in event.Digi_MuFilterHits :
        if mufilter_hit.GetSystem() == 3 :
            n_ds_hits += 1
            
            print(mufilter_hit.GetDetectorID())
            print(mufilter_hit.GetDetectorID()//1000)
            print(mufilter_hit.GetDetectorID()//1000-mufilter_hit.GetDetectorID()//10000*10)
            if mufilter_hit.GetDetectorID()//1000-mufilter_hit.GetDetectorID()//10000*10 >= 2 :
                n_ds_hits_last += 1

            if mufilter_hit.GetDetectorID() not in hit_list :
                n_ds_hits_unassigned += 1
    plot_vars["n_unassigned_hits"].append(n_ds_hits_unassigned)
    plot_vars["n_ds_hits"].append(n_ds_hits)
    plot_vars["n_last_ds_hits"].append(n_ds_hits_last)

    if plot_vars["n_true_mu"][-1] != plot_vars["n_reco_mu_raw"][-1] :
        print("Did not correctly identify muons. Expected {0}, reconstructed {1}".format(plot_vars["n_true_mu"][-1], plot_vars["n_reco_mu_raw"][-1]))
        print("File {0} event {1}".format(event.GetFile().GetPath(), event.GetReadEntry()-event.GetTreeOffset()[event.GetTreeNumber()]))
        misid_filename.append(event.GetFile().GetPath())
        misid_evtno.append(event.GetReadEntry()-event.GetTreeOffset()[event.GetTreeNumber()])
        misid_ntrue.append(plot_vars["n_true_mu"][-1])
        misid_nreco.append(plot_vars["n_reco_mu_raw"][-1])

    # At least one true muon in the event
    if len(true_muons) == 0 :
        continue

    cutTracker(signal_definition_true_cuts, "At least one true muon", weight)
    cutTracker(event_selection_reco_cuts, "At least one true muon", weight)
    
    if n_raw_reco_mu == 0 :
        continue
        
    cutTracker(event_selection_reco_cuts, "At least one reco muon", weight)

# Print some numbers
for k in plot_vars.keys() :
    plot_vars[k] = np.array(plot_vars[k])

for n_mu in range(0, 5) :
    print("True {0} muon events {1} correctly reconstructed {2}".format(n_mu, np.sum(plot_vars["w"][plot_vars["n_true_mu"] == n_mu]), np.sum(plot_vars["w"][np.logical_and(plot_vars["n_true_mu"] == n_mu, plot_vars["n_reco_mu_raw"] == n_mu)])))
print("Total events {0}".format(np.sum(plot_vars["w"])))

# Save misid info to text file
with open("misid_"+sample_name+".txt", "w") as f :
    for i in range(len(misid_filename)) :
        f.write("True {0} Reco {1} Filename {2} Event number {3}\n".format(misid_ntrue[i], misid_nreco[i], misid_filename[i], misid_evtno[i]))

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
plt.hist2d(plot_vars["n_true_mu"], plot_vars["n_reco_mu_raw"], weights = plot_vars["w"], range = ((0, 5), (0, 5)), bins = (5, 5))
plt.xlabel("Number of true muons")
plt.ylabel("Number of raw reconstructed muons")
plt.tight_layout()
plt.savefig("muon_multiplicity_raw_"+sample_name+".png")

plt.figure()
plt.hist2d(plot_vars["n_true_mu"], plot_vars["n_reco_mu_pruned"], weights = plot_vars["w"], range = ((0, 5), (0, 5)), bins = (5, 5))
plt.xlabel("Number of true muons")
plt.ylabel("Number of reconstructed muons")
plt.tight_layout()
plt.savefig("muon_multiplicity_pruned_"+sample_name+".png")

good_reco_mask = np.equal(plot_vars["n_true_mu"], plot_vars["n_reco_mu_raw"])

plt.figure()
plt.hist(plot_vars["n_unassigned_hits"][good_reco_mask], bins = 100, range = (0, 100), label = "Correctly reconstructed muon(s)", weights = plot_vars["w"][good_reco_mask], histtype = "step")
plt.hist(plot_vars["n_unassigned_hits"][~good_reco_mask], bins = 100, range = (0, 100), label = "Misidentified muon(s)", weights = plot_vars["w"][~good_reco_mask], histtype = "step")
plt.yscale("log")
plt.xlabel("Number of unassigned ds hits")

plt.figure()
plt.hist(plot_vars["n_ds_hits"][good_reco_mask], bins = 100, range = (0, 100), label = "Correctly reconstructed muon(s)", weights = plot_vars["w"][good_reco_mask], histtype = "step")
plt.hist(plot_vars["n_ds_hits"][~good_reco_mask], bins = 100, range = (0, 100), label = "Misidentified muon(s)", weights = plot_vars["w"][~good_reco_mask], histtype = "step")
plt.yscale("log")
plt.xlabel("Number of ds hits")

plt.figure()
plt.hist(plot_vars["n_last_ds_hits"][good_reco_mask], bins = 100, range = (0, 100), label = "Correctly reconstructed muon(s)", weights = plot_vars["w"][good_reco_mask], histtype = "step")
plt.hist(plot_vars["n_last_ds_hits"][~good_reco_mask], bins = 100, range = (0, 100), label = "Misidentified muon(s)", weights = plot_vars["w"][~good_reco_mask], histtype = "step")
plt.yscale("log")
plt.xlabel("Number of ds hits in layers 3 and 4")


plt.figure()
plt.hist(plot_vars["transverse_distance_reco_mu_shower_center"][good_reco_mask], bins = 75, range = (0, 75), label = "Correctly reconstructed muon(s)", weights = plot_vars["w"][good_reco_mask], histtype = "step")
plt.hist(plot_vars["transverse_distance_reco_mu_shower_center"][~good_reco_mask], bins = 75, range = (0, 75), label = "Misidentified muon(s)", weights = plot_vars["w"][~good_reco_mask], histtype = "step")
plt.yscale("log")
plt.show()

