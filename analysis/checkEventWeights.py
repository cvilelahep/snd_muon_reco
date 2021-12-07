import ROOT

import numpy as np

import matplotlib.pyplot as plt

ch = ROOT.TChain("cbmsim")
for i in range(100) :
    if i == 80 :
        continue
    ch.Add("/eos/home-c/cvilela/SND/neutrino/test_1k/{0}/sndLHC.Genie-TGeant4_dig.root".format(i))

n_entries = ch.GetEntries()

print("Got {0} entries".format(n_entries))

vtx = np.zeros((3, n_entries))
weight = np.zeros(n_entries)

for i, event in enumerate(ch) :
    if i % 100 == 0 :
        print("Reading entry {0}".format(i))
    for track in event.MCTrack :
        if track.GetMotherId() == -1 :
            vtx[0,i] = track.GetStartX()
            vtx[1,i] = track.GetStartY()
            vtx[2,i] = track.GetStartZ()
            
            weight[i] = track.GetWeight()


bins = 50
range_x = (-100, 100)
range_y = (-150, 70)
range_z = (-50, 50)

plt.figure()
plt.hist(vtx[0], label  = "Unweighted", bins = bins, range = range_x, histtype = "step")
plt.hist(vtx[0], label  = "Weighted", weights = weight, bins = bins, range = range_x, histtype = "step")
plt.xlabel("Vertex x [cm]")
plt.legend()

plt.figure()
plt.hist(vtx[1], label  = "Unweighted", bins = bins, range = range_y, histtype = "step")
plt.hist(vtx[1], label  = "Weighted", weights = weight, bins = bins, range = range_y, histtype = "step")
plt.legend()
plt.xlabel("Vertex y [cm]")

plt.figure()
plt.hist(vtx[2], label  = "Unweighted", histtype = "step", bins = bins, range = range_z)
plt.hist(vtx[2], label  = "Weighted", weights = weight, bins = bins, range = range_z, histtype = "step")
plt.legend()
plt.xlabel("Vertex Z [cm]")

plt.figure()
plt.subplot(2, 1, 1)
plt.hist2d(vtx[2], vtx[0],  bins = (bins, bins), range = (range_z, range_x))
plt.ylabel("Vertex x [cm]")
plt.title = "Unweighted"
plt.subplot(2, 1, 2)
plt.hist2d(vtx[2], vtx[0],  bins = (bins, bins), range = (range_z, range_x), weights = weight)
plt.ylabel("Vertex x [cm]")
plt.xlabel("Vertex z [cm]")
plt.title = "Weighted"

plt.figure()
plt.subplot(2, 1, 1)
plt.hist2d(vtx[2], vtx[1],  bins = (bins, bins), range = (range_z, range_y))
plt.ylabel("Vertex y [cm]")
plt.title = "Unweighted"
plt.subplot(2, 1, 2)
plt.hist2d(vtx[2], vtx[1],  bins = (bins, bins), range = (range_z, range_y), weights = weight)
plt.ylabel("Vertex y [cm]")
plt.xlabel("Vertex z [cm]")
plt.title = "Weighted"

plt.figure()
plt.subplot(2, 1, 1)
plt.hist2d(vtx[0], vtx[1],  bins = (bins, bins), range = (range_x, range_y))
plt.ylabel("Vertex y [cm]")
plt.title = "Unweighted"
plt.subplot(2, 1, 2)
plt.hist2d(vtx[0], vtx[1],  bins = (bins, bins), range = (range_x, range_y), weights = weight)
plt.ylabel("Vertex y [cm]")
plt.xlabel("Vertex x [cm]")
plt.title = "Weighted"

plt.show()
