import os.path

import numpy as np

import ROOT

import matplotlib.pyplot as plt

t_numu = ROOT.TChain("cbmsim")

numu_file_names = "/eos/home-c/cvilela/SND/neutrino/test_1k/{0}/sndLHC.Genie-TGeant4_dig.root"

n_muons = []
n_muons_thr = []
frac_pion_parent_thr = []
frac_kaon_parent_thr = []
frac_gamma_parent_thr = []
frac_nu_parent_thr = []
hard_muon_thr = 2.
muon_energies = []
muon_is_primary = []
muon_parent = []
nu_energies = []

for i in range(20) :
#for i in range(10) :
#for i in range(100) :
    if i == 80 :
        # Broken numu file
        continue
    print("Adding neutrino file {0}".format(i))
    if os.path.isfile(numu_file_names.format(i)) :
        t_numu.Add(numu_file_names.format(i))
    else :
        print("Skipping neutrino file {0}. Doesn't exist".format(i))

print("Done adding files")
print("Numu events", t_numu.GetEntries())

for i_event, event in enumerate(t_numu) :
    
    if event.MCTrack[0].GetStartX() < -47 :
        continue
    if event.MCTrack[0].GetStartX() > -8 :
        continue
    if event.MCTrack[0].GetStartY() < 15.5 :
        continue
    if event.MCTrack[0].GetStartY() > 54.5 :
        continue

    nu_energies.append( event.MCTrack[0].GetEnergy() )

    if nu_energies[-1] < 100 : 
        continue

    print(i_event)
    this_n_muons = 0
    this_n_muons_thr = 0
    this_frac_pion_thr = 0.
    this_frac_kaon_thr = 0.
    this_frac_gamma_thr = 0.
    this_frac_nu_thr = 0.
    for track in event.MCTrack :

        if track.GetStartX() < -47 :
            continue
        if track.GetStartX() > -8 :
            continue
        if track.GetStartY() < 15.5 :
            continue
        if track.GetStartY() > 54.5 :
            continue

        if abs(track.GetPdgCode()) != 13 :
            continue
        this_n_muons += 1
        muon_energies.append(track.GetEnergy())
        muon_is_primary.append(track.GetMotherId() == 0)
        muon_parent.append(event.MCTrack[track.GetMotherId()].GetPdgCode())
        if muon_energies[-1] > hard_muon_thr :
            this_n_muons_thr += 1
            if abs(muon_parent[-1]) == 14 :
                this_frac_nu_thr += 1
            elif abs(muon_parent[-1]) == 211 :
                this_frac_pion_thr += 1
            elif abs(muon_parent[-1]) == 321 :
                this_frac_kaon_thr += 1
            elif abs(muon_parent[-1]) == 22 :
                this_frac_gamma_thr += 1
    n_muons.append(this_n_muons)
    n_muons_thr.append(this_n_muons_thr)
    if this_n_muons_thr > 0 :
        frac_pion_parent_thr.append(this_frac_pion_thr/this_n_muons_thr)
        frac_kaon_parent_thr.append(this_frac_kaon_thr/this_n_muons_thr)
        frac_gamma_parent_thr.append(this_frac_gamma_thr/this_n_muons_thr)
        frac_nu_parent_thr.append(this_frac_nu_thr/this_n_muons_thr)
    else :
        frac_pion_parent_thr.append(0.)
        frac_kaon_parent_thr.append(0.)
        frac_gamma_parent_thr.append(0.)
        frac_nu_parent_thr.append(0.)

muon_energies = np.array(muon_energies)
muon_is_primary = np.array(muon_is_primary)
muon_parent = np.array(muon_parent)

frac_pion_parent_thr = np.array(frac_pion_parent_thr)
frac_kaon_parent_thr = np.array(frac_kaon_parent_thr)
frac_gamma_parent_thr = np.array(frac_gamma_parent_thr)
frac_nu_parent_thr = np.array(frac_nu_parent_thr)

n_muons_thr = np.array(n_muons_thr)

plt.figure()
plt.hist(n_muons, bins = 15, range = (0, 15))
plt.yscale('log')
plt.xlabel("Number of primary and secondary muons")
plt.tight_layout()
plt.savefig("multi_muon_multiplicity.png")

plt.figure()
plt.hist(nu_energies, bins = 100, range = (0, 5000))
plt.yscale('log')
plt.xlabel("Incoming neutrino energy [GeV]")
plt.tight_layout()
plt.savefig("multi_muon_nu_energy.png")


print(n_muons_thr.shape)
print(frac_pion_parent_thr.shape)

plt.figure()
plt.hist(n_muons_thr, bins = 5, range = (0, 5), histtype = 'step', label = "All")
plt.hist(n_muons_thr, bins = 5, range = (0, 5), histtype = 'step', weights = frac_gamma_parent_thr, label = "Gamma parent fraction")
plt.hist(n_muons_thr, bins = 5, range = (0, 5), histtype = 'step', weights = frac_pion_parent_thr, label = "Pion parent fraction")
plt.hist(n_muons_thr, bins = 5, range = (0, 5), histtype = 'step', weights = frac_kaon_parent_thr, label = "Kaon parent fraction")
plt.hist(n_muons_thr, bins = 5, range = (0, 5), histtype = 'step', weights = frac_nu_parent_thr, label = "Neutrino parent fraction")
plt.yscale('log')
plt.xlabel("Number of primary and secondary muons above 2 GeV")
plt.legend()
plt.tight_layout()
plt.savefig("multi_muon_multiplicity_thr.png")

plt.figure()
plt.hist(muon_energies, bins = 100, range = (0, 5000), histtype = 'step', label = "All")
plt.hist(muon_energies[muon_is_primary], bins = 100, range = (0, 5000), histtype = 'step', label = "Primary")
plt.hist(muon_energies[~muon_is_primary], bins = 100, range = (0, 5000), histtype = 'step', label = "Secondary")
plt.yscale('log')
plt.xlabel("Energy of primary and secondary muons [GeV]")
plt.legend()
plt.tight_layout()
plt.savefig("multi_muon_energy.png")

plt.figure()
plt.hist(muon_energies, bins = 500, range = (0, 5.), histtype = 'step', label = 'All')
plt.hist(muon_energies[muon_is_primary], bins = 500, range = (0, 5.), histtype = 'step', label = 'Primary')
plt.hist(muon_energies[~muon_is_primary], bins = 500, range = (0, 5.), histtype = 'step', label = 'Secondary')
plt.yscale('log')
plt.xlabel("Energy of primary and secondary muons [GeV]")
plt.legend()
plt.tight_layout()
plt.savefig("multi_muon_energy_zoomed.png")

plt.show()
