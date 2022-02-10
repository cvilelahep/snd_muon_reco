import ROOT
import numpy as np
import glob

from collections import OrderedDict

import matplotlib.pyplot as plt

# To store detector positions
a = ROOT.TVector3()
b = ROOT.TVector3()

# Target dims for fiducial volume. Approximately estimated from event display
x_extent = [-47, -6]
y_extent = [18, 53]
z_extent = [288, 354]

#x_extent = [-47, -8]
#y_extent = [15.5, 54.5]
#z_extent = [-25, 25]

def fiducialVolume(x, y, z) :
    # Target station dimensions. Consider making this more precise.
    print(x, y, z)
    if x < x_extent[0] :
        return False
    if x > x_extent[1] :
        return False
    if y < y_extent[0] :
        return False
    if y > y_extent[1] :
        return False
    if z < z_extent[0] :
        return False
    if z > z_extent[1] :
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
    
    return (fiducialVolume(track.GetStartX(), track.GetStartY(), track.GetStartZ()), weight, track.GetStartX(), track.GetStartY(), track.GetStartZ())


input_filename_expr = ["/eos/home-c/cvilela/SND_FEB_7/test/sndLHC.Genie-TGeant4"]
# Sort out input files
input_digi_filenames = []
input_reco_filenames = []
for expr in input_filename_expr :
    for fname in glob.glob(expr+"_dig.root") :
        input_digi_filenames.append(fname)
#    for fname in glob.glob(expr+"_dig_muonReco.root") :
    for fname in glob.glob(expr+"_dig_muonReco_tol0.root") :
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

w = []

plot_vars = {}
plot_vars["vtx_x"] = []
plot_vars["vtx_y"] = []
plot_vars["vtx_z"] = []

plot_vars["true_angle_zx"] = []
plot_vars["true_angle_zy"] = []

plot_vars["hough_angle_zx"] = []
plot_vars["hough_angle_zy"] = []

plot_vars["kalman_angle_zx"] = []
plot_vars["kalman_angle_zy"] = []

n_true_1mu = 0.
n_reco_1mu = 0.

for i_event, event in enumerate(tree_reco) :
    print(i_event)
    inFV, weight, x, y, z = getNuWeightFV(event)

    if not inFV :
        print("Not in FV")
        continue

    true_muons = trueMuons(event)
    if len(true_muons) != 1 :
        continue 

    n_true_1mu += weight

    hough_muons = event.Reco_MuonTracks
    
    if len(hough_muons) != 1 :
        continue
    n_reco_1mu += weight

    hough_muon = event.Reco_MuonTracks[0]

    w.append(weight)
    plot_vars["vtx_x"].append(x)
    plot_vars["vtx_y"].append(y)
    plot_vars["vtx_z"].append(z)


    true_muon = event.MCTrack[true_muons[0]]

    plot_vars["true_angle_zx"].append(np.arctan2( true_muon.GetPx(),
                                                  true_muon.GetPz() ))
    plot_vars["true_angle_zy"].append(np.arctan2( true_muon.GetPy(),
                                                  true_muon.GetPz() ))

    plot_vars["hough_angle_zx"].append(np.arctan2( (hough_muon.getStop().X() - hough_muon.getStart().X()),
                                                   (hough_muon.getStop().Z() - hough_muon.getStart().Z()) ))
    plot_vars["hough_angle_zy"].append(np.arctan2( (hough_muon.getStop().Y() - hough_muon.getStart().Y()),
                                                   (hough_muon.getStop().Z() - hough_muon.getStart().Z()) ))
    print("Got Hough")
#    for kalman_muon in event.Reco_KalmanTracks :
        
    kalman_muon = event.Reco_KalmanTracks[0]
#        print(kalman_muon)
#    print(dir(kalman_muon))
#        print(kalman_muon.getNumPointsWithMeasurement())
#        print(kalman_muon.getFittedState())
    plot_vars["kalman_angle_zx"].append(np.arctan2( kalman_muon.getFittedState().getMom().x(),
                                                    kalman_muon.getFittedState().getMom().z() ))
    plot_vars["kalman_angle_zy"].append(np.arctan2( kalman_muon.getFittedState().getMom().y(),
                                                    kalman_muon.getFittedState().getMom().z() ))
    print("Got Kalman")
    
    
w = np.array(w)
for k in plot_vars.keys() :
    plot_vars[k] = np.array(plot_vars[k])
                                    


print("1mu efficiency {0}".format(n_reco_1mu/n_true_1mu))
print("MAKING PLOTS!")
plt.figure()
plt.hist(plot_vars["vtx_x"], weights = w, bins = 50, histtype = 'step', density = True)
plt.hist(plot_vars["vtx_x"], bins = 50, histtype = 'step', density = True)
plt.xlabel("True vertex X [cm]")
#plt.legend()
plt.tight_layout()

plt.figure()
plt.hist(plot_vars["vtx_y"], weights = w, bins = 50, histtype = 'step', density = True)
plt.hist(plot_vars["vtx_y"], bins = 50, histtype = 'step', density = True)
plt.xlabel("True vertex Y [cm]")
#plt.legend()
plt.tight_layout()

plt.figure()
plt.hist(plot_vars["vtx_z"], weights = w, bins = 50, histtype = 'step', density = True)
plt.hist(plot_vars["vtx_z"], bins = 50, histtype = 'step', density = True)
plt.xlabel("True vertex Z [cm]")
#plt.legend()
plt.tight_layout()

plt.figure()
plt.hist((plot_vars["hough_angle_zx"]-plot_vars["true_angle_zx"]), weights = w, bins = 100, histtype = "step", label = "Hough")
plt.hist((plot_vars["kalman_angle_zx"]-plot_vars["true_angle_zx"]), weights = w, bins = 100, histtype = "step", label = "Kalman")
plt.xlabel("ZX reco - true angle [rad]")
plt.legend()
plt.tight_layout()

plt.figure()
plt.hist((plot_vars["hough_angle_zy"]-plot_vars["true_angle_zy"]), weights = w, bins = 100, histtype = "step", label = "Hough")
plt.hist((plot_vars["kalman_angle_zy"]-plot_vars["true_angle_zy"]), weights = w, bins = 100, histtype = "step", label = "Kalman")
plt.xlabel("ZY reco - true angle [rad]")
plt.legend()
plt.tight_layout()


plt.show()

