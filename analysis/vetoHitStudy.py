import os.path

import numpy as np

import ROOT

# Geometry stuff
from rootpyPickler import Unpickler

fgeo = ROOT.TFile.Open("/eos/home-c/cvilela/SND/muon/muons_VCdown_IR1-LHC/10/geofile_full.Ntuple-TGeant4.root")
upkl    = Unpickler(fgeo)
snd_geo = upkl.load('ShipGeo')
run = ROOT.FairRunSim()

import shipLHC_conf as sndDet_conf
modules = sndDet_conf.configure(run,snd_geo)
sGeo = fgeo.FAIRGeom
modules['Scifi'].SiPMmapping()

background_file_names = "/eos/home-c/cvilela/SND/muon/muons_VCdown_IR1-LHC/{0}/sndLHC.Ntuple-TGeant4_dig.root"
numu_file_names = "/eos/home-c/cvilela/SND/neutrino/test_1k/{0}/sndLHC.Genie-TGeant4_dig.root"

t_background = ROOT.TChain("cbmsim")
t_numu = ROOT.TChain("cbmsim")

nu_veto_hit_multiplicity = []
nu_veto_hit_multiplicity_thr = []
nu_veto_hit_total_energy = []
nu_veto_hit_energy_per_hit = []
nu_total_energy = []
nu_tof = []
nu_tof_true = []
nu_weight = []
nu_z = []

background_veto_hit_multiplicity = []
background_veto_hit_multiplicity_thr = []
background_veto_hit_total_energy = []
background_veto_hit_energy_per_hit = []
background_total_energy = []
background_tof = []

background_tof_true = []

tof_cut = 0.5
veto_hit_threshold = 0.0015
minimum_energy_cut = 0.015

#for i in range(10) :
for i in range(100) :
    if i == 80 :
        # Broken numu file
        continue
    if os.path.isfile(background_file_names.format(i)) :
        print("Adding backfround file {0}".format(i))
        t_background.Add(background_file_names.format(i))
    else :
        print("Skipping background file {0}. Doesn't exist".format(i))
    print("Adding neutrino file {0}".format(i))
    if os.path.isfile(numu_file_names.format(i)) :
        t_numu.Add(numu_file_names.format(i))
    else :
        print("Skipping neutrino file {0}. Doesn't exist".format(i))

print("Done adding files")
print("Numu events", t_numu.GetEntries())
print("Background events", t_background.GetEntries())

a = ROOT.TVector3()
b = ROOT.TVector3()

for i, nu_event in enumerate(t_numu) :
    if i%10 == 0 :
        print("Running nu event {0}".format(i))

    try :
        if len(nu_event.Digi_ScifiHits) == 0 :
            continue
    except AttributeError :
        continue

    for track in nu_event.MCTrack :
        if track.GetMotherId() == -1 :
            break

    weight = track.GetWeight()
    if track.GetStartX() < -47 :
        continue
    if track.GetStartX() > -8 :
        continue
    if track.GetStartY() < 15.5 :
        continue
    if track.GetStartY() > 54.5 :
        continue


    n_hits = 0
    n_hits_thr = 0
    tot_veto_energy = 0
    tot_energy = 0

    earliest_scifi_hit_t = 1e10
    earliest_scifi_hit_t_true = 1e10
    earliest_scifi_hit_z = -500

    for point in nu_event.ScifiPoint :
        t = point.GetTime()
        t_smear = np.random.normal(t, 0.25)

        if t_smear < earliest_scifi_hit_t :
            earliest_scifi_hit_t = t_smear
            modules["Scifi"].GetSiPMPosition(point.GetDetectorID(), a, b)
            earliest_scifi_hit_z = a.Z()
        if t < earliest_scifi_hit_t_true :
            earliest_scifi_hit_t_true = t
    
    earliest_veto_hit_t = 1e10
    earliest_veto_hit_t_true = 1e10
    earliest_veto_hit_z = -500

    for point in nu_event.MuFilterPoint :

        if np.floor(point.GetDetectorID()/10000) != 1 :
            continue

        t = point.GetTime()
        t_smear = np.random.normal(t, 0.1)

        if t_smear < earliest_veto_hit_t :
            earliest_veto_hit_t = t_smear
            modules["MuFilter"].GetPosition(point.GetDetectorID(), a, b)
            earliest_veto_hit_z = a.Z()
        if t < earliest_veto_hit_t_true :
            earliest_veto_hit_t_true = t


    for hit in nu_event.Digi_MuFilterHits :
        tot_energy += hit.GetEnergy()
        if hit.GetSystem() == 1 :
            n_hits += 1
            tot_veto_energy += hit.GetEnergy()
            if hit.GetEnergy() > veto_hit_threshold :
                n_hits_thr += 1

#    for hit in nu_event.Digi_ScifiHits :
#        tot_energy += hit.GetEnergy()

    if tot_energy < minimum_energy_cut  :
        continue

    nu_weight.append(weight)    
    nu_z.append(track.GetStartZ())

    nu_tof.append(earliest_veto_hit_t - earliest_scifi_hit_t)
    nu_tof_true.append(earliest_veto_hit_t_true - earliest_scifi_hit_t_true)


    nu_veto_hit_multiplicity.append(n_hits)
    nu_veto_hit_multiplicity_thr.append(n_hits_thr)
    nu_veto_hit_total_energy.append(tot_veto_energy)
    nu_total_energy.append(tot_energy)
    if n_hits > 0 :
        nu_veto_hit_energy_per_hit.append(tot_veto_energy/n_hits)
    else :
        nu_veto_hit_energy_per_hit.append(-1)

print("Done numus")

for i, nu_event in enumerate(t_background) :

    if i%10 == 0 :
        print("Running mu event {0}".format(i))
    
    if len(nu_event.Digi_ScifiHits) == 0 :
        continue

    n_hits = 0
    n_hits_thr = 0
    tot_veto_energy = 0
    tot_energy = 0
    for hit in nu_event.Digi_MuFilterHits :
        tot_energy += hit.GetEnergy()
        if hit.GetSystem() == 1 :
            n_hits += 1
            tot_veto_energy += hit.GetEnergy()
            if hit.GetEnergy() > veto_hit_threshold :
                n_hits_thr += 1

    earliest_scifi_hit_t = 1e10
    earliest_scifi_hit_t_true = 1e10
    earliest_scifi_hit_z = -500

    for point in nu_event.ScifiPoint :
        t = point.GetTime()
        t_smear = np.random.normal(t, 0.25)

        if t_smear < earliest_scifi_hit_t :
            earliest_scifi_hit_t = t_smear
            modules["Scifi"].GetSiPMPosition(point.GetDetectorID(), a, b)
            earliest_scifi_hit_z = a.Z()
        if t < earliest_scifi_hit_t_true :
            earliest_scifi_hit_t_true = t
    
    earliest_veto_hit_t = 1e10
    earliest_veto_hit_t_true = 1e10
    earliest_veto_hit_z = -500

    for point in nu_event.MuFilterPoint :

        if np.floor(point.GetDetectorID()/10000) != 1 :
            continue

        t = point.GetTime()
        t_smear = np.random.normal(t, 0.1)

        if t_smear < earliest_veto_hit_t :
            earliest_veto_hit_t = t_smear
            modules["MuFilter"].GetPosition(point.GetDetectorID(), a, b)
            earliest_veto_hit_z = a.Z()
        if t < earliest_veto_hit_t_true :
            earliest_veto_hit_t_true = t

    if tot_energy < minimum_energy_cut  :
        continue

    background_tof.append(earliest_veto_hit_t - earliest_scifi_hit_t)
    background_tof_true.append(earliest_veto_hit_t_true - earliest_scifi_hit_t_true)

#    for hit in nu_event.Digi_ScifiHits :
#        tot_energy += hit.GetEnergy()

    background_veto_hit_multiplicity.append(n_hits)
    background_veto_hit_multiplicity_thr.append(n_hits_thr)
    background_veto_hit_total_energy.append(tot_veto_energy)
    background_total_energy.append(tot_energy)
    if n_hits > 0 :
        background_veto_hit_energy_per_hit.append(tot_veto_energy/n_hits)

print("Done muons")

nu_veto_hit_multiplicity = np.array(nu_veto_hit_multiplicity)
nu_veto_hit_multiplicity_thr = np.array(nu_veto_hit_multiplicity_thr)
nu_veto_hit_total_energy = np.array(nu_veto_hit_total_energy)
nu_total_energy = np.array(nu_total_energy)
nu_veto_hit_energy_per_hit = np.array(nu_veto_hit_energy_per_hit)
nu_tof = np.array(nu_tof)
nu_tof_true = np.array(nu_tof_true)
nu_weight = np.array(nu_weight)
nu_z = np.array(nu_z)
background_veto_hit_multiplicity = np.array(background_veto_hit_multiplicity)
background_veto_hit_multiplicity_thr = np.array(background_veto_hit_multiplicity_thr)
background_veto_hit_total_energy = np.array(background_veto_hit_total_energy)
background_total_energy = np.array(background_total_energy)
background_veto_hit_energy_per_hit = np.array(background_veto_hit_energy_per_hit)                    
background_tof = np.array(background_tof)
background_tof_true = np.array(background_tof_true)


import matplotlib.pyplot as plt
plt.figure()
#plt.hist(nu_veto_hit_multiplicity, histtype = "step", bins = 20, range = (0, 20),  label = r"$\nu_{\mu}$CCDIS")
plt.hist(nu_veto_hit_multiplicity, histtype = "step", bins = 20, range = (0, 20),  label = r"$\nu_{\mu}$CCDIS", weights = nu_weight/np.sum(nu_weight))
plt.hist(background_veto_hit_multiplicity, histtype = "step", bins = 20, range = (0, 20),  label = "Background muons", weights = [1./len(background_veto_hit_multiplicity)]*len(background_veto_hit_multiplicity))
plt.legend()
plt.xlabel("Veto hit multiplicity")
plt.tight_layout()
plt.savefig("snd_veto_multiplicity.png")

plt.figure()
#plt.hist(nu_veto_hit_multiplicity, histtype = "step", bins = 20, range = (0, 20),  label = r"$\nu_{\mu}$CCDIS")
plt.hist(nu_veto_hit_multiplicity_thr, histtype = "step", bins = 20, range = (0, 20),  label = r"$\nu_{\mu}$CCDIS", weights = nu_weight/np.sum(nu_weight))
plt.hist(background_veto_hit_multiplicity_thr, histtype = "step", bins = 20, range = (0, 20),  label = "Background muons", weights = [1./len(background_veto_hit_multiplicity_thr)]*len(background_veto_hit_multiplicity_thr))
plt.legend()
plt.xlabel("Veto hit multiplicity (1.0 MeV threshold)")
plt.tight_layout()
plt.savefig("snd_veto_multiplicity_thr.png")

plt.figure()
plt.hist(nu_veto_hit_multiplicity, histtype = "step", bins = 20, range = (0, 20),  label = r"$\nu_{\mu}$CCDIS", weights = nu_weight/np.sum(nu_weight), linestyle = '--')
plt.hist(background_veto_hit_multiplicity, histtype = "step", bins = 20, range = (0, 20),  label = "Background muons", weights = [1./len(background_veto_hit_multiplicity)]*len(background_veto_hit_multiplicity), linestyle = '--')
plt.hist(nu_veto_hit_multiplicity_thr, histtype = "step", bins = 20, range = (0, 20), weights = nu_weight/np.sum(nu_weight), color = 'tab:blue', label = r"$\nu_{\mu}$CCDIS (1.0 MeV threshold)")
plt.hist(background_veto_hit_multiplicity_thr, histtype = "step", bins = 20, range = (0, 20), weights = [1./len(background_veto_hit_multiplicity_thr)]*len(background_veto_hit_multiplicity_thr), color = 'tab:orange', label = "Background muons (1.0 MeV threshold)")
plt.xlabel("Veto hit multiplicity")
plt.legend()
plt.tight_layout()
plt.savefig("snd_veto_multiplicity_overlay.png")


plt.figure()
#plt.hist(nu_veto_hit_total_energy, histtype = "step", bins = 50, range = (0, 0.02),  label = r"$\nu_{\mu}$CCDIS")
plt.hist(nu_veto_hit_total_energy, histtype = "step", bins = 50, range = (0, 0.02),  label = r"$\nu_{\mu}$CCDIS", weights = nu_weight/np.sum(nu_weight))
plt.hist(background_veto_hit_total_energy, histtype = "step", bins = 50, range = (0, 0.02),  label = "Background muons", weights = [1./len(background_veto_hit_total_energy)]*len(background_veto_hit_total_energy))
plt.legend()
plt.xlabel("Total energy deposited in veto system [GeV]")
plt.tight_layout()
plt.savefig("snd_veto_total_energy.png")

plt.figure()
#plt.hist(nu_veto_hit_total_energy, histtype = "step", bins = 50, range = (0, 0.02),  label = r"$\nu_{\mu}$CCDIS")
plt.hist(nu_total_energy, histtype = "step", bins = 100, range = (0, 0.5),  label = r"$\nu_{\mu}$CCDIS", weights = nu_weight/np.sum(nu_weight))
plt.hist(background_total_energy, histtype = "step", bins = 100, range = (0, 0.5),  label = "Background muons", weights = [1./len(background_total_energy)]*len(background_total_energy))
plt.legend()
plt.xlabel("Total energy deposited in muon and veto systems [GeV]")
plt.tight_layout()
plt.savefig("snd_total_energy.png")
plt.figure()

#plt.hist(nu_veto_hit_total_energy, histtype = "step", bins = 50, range = (0, 0.02),  label = r"$\nu_{\mu}$CCDIS")
plt.subplot(1, 2, 1)
plt.hist2d(nu_total_energy, nu_veto_hit_multiplicity, bins = (50, 5), range = ((0, .25), (0, 5)), weights = nu_weight/np.sum(nu_weight))
plt.subplot(1, 2, 2)
plt.hist2d(background_total_energy, background_veto_hit_multiplicity, bins = (50, 5), range = ((0, .25), (0, 5)), weights = [1./len(background_total_energy)]*len(background_total_energy))
#plt.xlabel("Total energy deposited in muon and veto systems [GeV]")
plt.tight_layout()
plt.savefig("snd_total_energy_veto_hits.png")


plt.figure()
#plt.hist(nu_veto_hit_total_energy, histtype = "step", bins = 50, range = (0, 0.02),  label = r"$\nu_{\mu}$CCDIS")
plt.hist(nu_veto_hit_energy_per_hit[nu_veto_hit_energy_per_hit>=0], histtype = "step", bins = 50, range = (0, 0.02),  label = r"$\nu_{\mu}$CCDIS", weights = nu_weight[nu_veto_hit_energy_per_hit>=0]/np.sum(nu_weight[nu_veto_hit_energy_per_hit>=0]))
plt.hist(background_veto_hit_energy_per_hit[background_veto_hit_energy_per_hit>=0], histtype = "step", bins = 50, range = (0, 0.02),  label = "Background muons", weights = [1./len(background_veto_hit_energy_per_hit[background_veto_hit_energy_per_hit>=0])]*len(background_veto_hit_energy_per_hit[background_veto_hit_energy_per_hit>=0]))
plt.legend()
plt.xlabel("Average energy per veto system hit [GeV]")
plt.tight_layout()
plt.savefig("snd_veto_energy_per_hit.png")

plt.figure()
plt.subplot(1, 2, 1)
plt.hist2d(nu_veto_hit_total_energy, nu_veto_hit_multiplicity, bins = (50, 20), range = ((0, 0.02), (0, 20)), weights = nu_weight/np.sum(nu_weight))
plt.subplot(1, 2, 2)
plt.hist2d(background_veto_hit_total_energy, background_veto_hit_multiplicity, bins = (50, 20), range = ((0, 0.02), (0, 20)))
plt.tight_layout()
plt.savefig("snd_veto_multiplicity_vs_total_veto_energy.png")

plt.figure()
plt.hist(nu_tof[nu_veto_hit_multiplicity_thr>0], histtype = "step", bins = 50, range = (-5, 5), label = r"$\nu_{\mu}$CCDIS", weights = nu_weight[nu_veto_hit_multiplicity_thr>0]/np.sum(nu_weight[nu_veto_hit_multiplicity_thr>0]))
plt.hist(background_tof[background_veto_hit_multiplicity_thr>0], histtype = "step", bins = 50, range = (-5, 5), label = "Background muons", weights = [1./len(background_tof[background_veto_hit_multiplicity_thr>0])]*len(background_tof[background_veto_hit_multiplicity_thr>0]))
plt.hist(nu_tof_true[nu_veto_hit_multiplicity_thr>0], histtype = "step", bins = 50, range = (-5, 5), label = r"$\nu_{\mu}$CCDIS (no t smearing)", weights = nu_weight[nu_veto_hit_multiplicity_thr>0]/np.sum(nu_weight[nu_veto_hit_multiplicity_thr>0]), color = "tab:blue", linestyle = "--")
plt.hist(background_tof_true[background_veto_hit_multiplicity_thr>0], histtype = "step", bins = 50, range = (-5, 5), label = "Background muons (no t smearing)", weights = [1./len(background_tof[background_veto_hit_multiplicity_thr>0])]*len(background_tof[background_veto_hit_multiplicity_thr>0]), color = "tab:orange", linestyle = "--")
plt.legend()
plt.xlabel("Earliest veto hit time minus earliest Scifi hit time [ns]")
plt.tight_layout()
plt.savefig("snd_veto_scifi_tof.png")

plt.figure()
n, bins, patches = plt.hist(nu_z, histtype = "step", bins = 50, range = (-40, 40), label = "All events", weights = nu_weight)
plt.hist(nu_z[nu_veto_hit_multiplicity_thr==0], histtype = "step", bins = 50, range = (-40, 40), label = "Events with no veto hits", weights = nu_weight[nu_veto_hit_multiplicity_thr==0])
plt.hist(nu_z[np.logical_and(nu_veto_hit_multiplicity_thr>0, nu_tof > tof_cut)], histtype = "step", bins = 50, range = (-40, 40), label = "Events with veto hits and internal TOF", weights = nu_weight[np.logical_and(nu_veto_hit_multiplicity_thr>0, nu_tof > tof_cut)])
plt.hist(nu_z[np.logical_or(nu_veto_hit_multiplicity_thr == 0, np.logical_and(nu_veto_hit_multiplicity_thr > 0, nu_tof > tof_cut))], histtype = "step", bins = 50, range = (-40, 40), label = "Events with no veto hits or with internal TOF", weights = nu_weight[np.logical_or(nu_veto_hit_multiplicity_thr==0, np.logical_and(nu_veto_hit_multiplicity_thr>0, nu_tof > tof_cut))])
plt.ylim(0, np.max(n)*1.5)
plt.legend()
plt.xlabel("True neutrino vertex Z [cm]")
plt.tight_layout()
plt.savefig("snd_veto_vtx_z.png")

print("NEUTRINOS")
print("TOTAL {0}, No veto hits above threshold {1}, At least one veto hit above threshold, but timing inconsistent with neutrino {2}".format(np.sum(nu_weight), np.sum(nu_weight[nu_veto_hit_multiplicity_thr==0]), np.sum(nu_weight[np.logical_and(nu_veto_hit_multiplicity_thr>0, nu_tof > tof_cut)])))
print("MUONS")
print("TOTAL {0}, No veto hits above threshold {1}, At least one veto hit above threshold, but timing inconsistent with neutrino {2}".format(len(background_veto_hit_multiplicity_thr), np.sum(background_veto_hit_multiplicity_thr==0), np.sum(np.logical_and(background_veto_hit_multiplicity_thr>0, background_tof > tof_cut))))

plt.show()
