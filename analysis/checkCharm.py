import ROOT
import glob
from collections import Counter

import matplotlib.pyplot as plt

import numpy as np
import scipy.optimize

input_filename_expr = ["/eos/user/c/cvilela/SND_ANALYSIS/neutrino/numuDefault/*/sndLHC.Genie-TGeant4_dig.root",
                       "/eos/user/c/cvilela/SND_ANALYSIS/neutrino/nueDefault/*/sndLHC.Genie-TGeant4_dig.root"]

input_digi_filenames = []
for expr in input_filename_expr :
    for fname in glob.glob(expr) :
        input_digi_filenames.append(fname)

events = ROOT.TChain("cbmsim")
for fname in input_digi_filenames :
    print("Adding {0}".format(fname))
    events.Add(fname)

counter = 0
counter_charm = 0
all_charm_tracks_PDG = []
charm_daughters = {}
charm_distances = {}
charm_mother_z = {}
charm_daughter_z = {}
charm_times = {}

for i_event, event in enumerate(events) :
    if i_event > 50000 :
        break
    counter += 1

    charm_tracks_ID = []
    charm_tracks_PDG = []
    charm_start_x = []
    charm_start_y = []
    charm_start_z = []
    charm_start_t = []
    for i_track, track in enumerate(event.MCTrack) :
        if abs(track.GetPdgCode()) in [411, 421, 4122, 431] :
            if track.GetMotherId() == 0 :
                if i_track not in charm_tracks_ID :
                    charm_tracks_ID.append(i_track)
                    charm_tracks_PDG.append(track.GetPdgCode())
                    charm_start_x.append(track.GetStartX())
                    charm_start_y.append(track.GetStartY())
                    charm_start_z.append(track.GetStartZ())
                    charm_start_t.append(track.GetStartT())
                
    if len(charm_tracks_ID) :
        print("Event {0}".format(i_event))
        for i_charm in range(len(charm_tracks_ID)) :
            print("Track ID {0} PDG Code {1}".format(charm_tracks_ID[i_charm], charm_tracks_PDG[i_charm]))
        all_charm_tracks_PDG += charm_tracks_PDG
        counter_charm += 1
    for i_track, track in enumerate(event.MCTrack) :
        if track.GetMotherId() in charm_tracks_ID :
            mother_pdg = charm_tracks_PDG[charm_tracks_ID.index(track.GetMotherId())]
            mother_x = charm_start_x[charm_tracks_ID.index(track.GetMotherId())]
            mother_y = charm_start_y[charm_tracks_ID.index(track.GetMotherId())]
            mother_z = charm_start_z[charm_tracks_ID.index(track.GetMotherId())]
            mother_t = charm_start_t[charm_tracks_ID.index(track.GetMotherId())]

            if mother_pdg not in charm_mother_z :
                charm_mother_z[mother_pdg] = [mother_z]
            else :
                charm_mother_z[mother_pdg].append(mother_z)
            
            if mother_pdg not in charm_daughters :
                charm_daughters[mother_pdg] = [track.GetPdgCode()]
                charm_distances[mother_pdg] = [ ( (mother_x - track.GetStartX())**2
                                                  + (mother_y - track.GetStartY())**2
                                                  + (mother_z - track.GetStartZ())**2)**0.5 ]
                charm_times[mother_pdg] = [track.GetStartT() - mother_t]

                charm_daughter_z[mother_pdg] = [track.GetStartZ()]
            else :
                charm_daughters[mother_pdg].append(track.GetPdgCode())
                charm_distances[mother_pdg].append( ( (mother_x - track.GetStartX())**2
                                                      + (mother_y - track.GetStartY())**2
                                                      + (mother_z - track.GetStartZ())**2)**0.5 )
                charm_times[mother_pdg].append(track.GetStartT() - mother_t)
                charm_daughter_z[mother_pdg].append(track.GetStartZ())

            print("T mother {0} daughter {1}".format(mother_t, track.GetStartT()))
            print("Z mother {0} daughter {1}".format(mother_z, track.GetStartZ()))

print("Fraction of events with charm production: {0}".format(counter_charm/float(counter)))
charm_counts = Counter(all_charm_tracks_PDG)
for pdg, n in charm_counts.items() :
    print("Fraction of events with {0}: {1}".format(pdg, n/float(counter)))

print("Branching ratios")

plt.figure()
plt.bar(range(len(charm_counts)), [float(n)/len(all_charm_tracks_PDG) for n in charm_counts.values()])
plt.xticks(range(len(charm_counts)), charm_counts.keys(), rotation = 90)
plt.xlabel("Primary charm hadron")
plt.tight_layout()
plt.savefig("charm_mother_PDG.png")

for mother_pdg, daughters in charm_daughters.items() :
    
    daughter_counts = Counter(daughters)
    print("Mother PDG {0}".format(mother_pdg))
    for pdg, n in daughter_counts.items() :
        print("Fraction of {0} daughters: {1}".format(pdg, float(n)/len(daughters)))
    plt.figure()
    plt.bar(range(len(daughter_counts)), [float(n)/len(daughters) for n in daughter_counts.values()])
    plt.xticks(range(len(daughter_counts)), daughter_counts.keys(), rotation = 90)
    plt.title("Fraction of {0} daughter particles".format(mother_pdg))
    plt.xlabel("Daughter PDG")
    plt.tight_layout()
    plt.savefig("charm_mother_daughter_PDG_"+str(mother_pdg)+".png")

def f_exp(x, a, b) :
    return a*np.exp(-x/b)

for mother_pdg, d in charm_distances.items() :
    plt.figure()
    n, bins, patches = plt.hist(d, bins = 100, range = (0, 5))
    x_centers = np.add(bins[:-1], bins[1:])/2.
    popt, pcov = scipy.optimize.curve_fit(f_exp, x_centers, n)
    plt.plot(x_centers, f_exp(x_centers, popt[0], popt[1]), label = "<d> = {0:.2g} cm".format(popt[1]))
    plt.legend()
    plt.xlabel("Distance travelled by charm hadron [cm]")
    plt.title("Charm hadron PDG: {0}".format(mother_pdg))
    plt.tight_layout()
    plt.savefig("charm_mother_daugther_distance_"+str(mother_pdg)+".png")
#    plt.yscale('log')

#for mother_pdg, t in charm_times.items() :
#    plt.figure()
#    plt.hist(t)
#    plt.xlabel("Charm hadron livetime")
#    plt.title("Charm hadron PDG: {0}".format(mother_pdg))
#
#for mother_pdg in charm_daughter_z.keys() :
#    plt.figure()
#    plt.hist(charm_mother_z[mother_pdg], histtype = "step", label = "Mother")
#    plt.hist(charm_daughter_z[mother_pdg], histtype = "step", label = "Daughter")
#    plt.xlabel("Starting position Z")
#    plt.title("Charm hadron PDG: {0}".format(mother_pdg))
    

plt.show()
