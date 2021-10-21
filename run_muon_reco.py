import argparse
import numpy as np
import pickle

import ROOT

from hough_transform import hough

# Geometry stuff
from rootpyPickler import Unpickler

parser = argparse.ArgumentParser()
parser.add_argument("input_file")
parser.add_argument("geo_file")
parser.add_argument("output_file")

args = parser.parse_args()

print("Input: ", args.input_file)
print("Output: ", args.output_file)
print("Geometry: ", args.geo_file)

fgeo = ROOT.TFile.Open(args.geo_file)
upkl    = Unpickler(fgeo)
snd_geo = upkl.load('ShipGeo')
run = ROOT.FairRunSim()

import shipLHC_conf as sndDet_conf
modules = sndDet_conf.configure(run,snd_geo)
sGeo = fgeo.FAIRGeom
modules['Scifi'].SiPMmapping()

events = ROOT.TChain("cbmsim")

events.Add(args.input_file)

# Reco parameters:
# Maximum absolute value of reconstructed angle
max_angle = np.pi*0.1
# Number of bins per Hough accumulator axis
n = 1000
# Number of random throws per hit
n_random = 5
# MuFilter weight. Muon filter hits are thrown more times than scifi
muon_weight = 100
# Minimum number of hits in each of the downstream muon filter views to try to reconstruct a muon
min_hits = 3

# Initialize Hough transforms for both views:
h_ZX = hough(n, [-80, 0], n, [-max_angle+np.pi/2., max_angle+np.pi/2.], False)
h_ZY = hough(n, [0, 80], n, [-max_angle+np.pi/2., max_angle+np.pi/2.], False)

# List to store output
out_data = []

# Vectors to read hit positions
a = ROOT.TVector3()
b = ROOT.TVector3()

# Loop over events
for i_event, event in enumerate(events) :

    this_data = []

#    if i_event % 100 == 0 :
    print("Running event {0}".format(i_event))

    # Read hits
    # For downstream muon filter hits
    mu_ds = {"pos" : [[], [], []], 
             "d" : [[], [], []], 
             "vert" : [], 
             "index" : [],
             "system" : []}

    # For upstream muon filter hits
    mu_us = {"pos" : [[], [], []], 
             "d" : [[], [], []], 
             "vert" : [],
             "index" : [],
             "system" : [] }

    # For scifi hits
    scifi = {"pos" : [[], [], []], 
             "d" : [[], [], []], 
             "vert" : [],
             "index" : [], 
             "system" : []}

    # Loop through hits
    for i_hit, muFilterHit in enumerate(event.Digi_MuFilterHits) :
        
        # Don't use veto for fitting
        if muFilterHit.GetSystem() == 1 :
            continue
        elif muFilterHit.GetSystem() == 2 :
            mu = mu_us
        elif muFilterHit.GetSystem() == 3 :
            mu = mu_ds
        else :
            print("WARNING! Unknown MuFilter system!!")

        modules['MuFilter'].GetPosition(muFilterHit.GetDetectorID(), a, b)
        mu["pos"][0].append(a.X())
        mu["pos"][1].append(a.Y())
        mu["pos"][2].append(a.Z())

        mu["vert"].append(muFilterHit.isVertical())
        mu["system"].append(muFilterHit.GetSystem())
        
        mu["d"][0].append(1.)
        mu["d"][2].append(1.)
        
        mu["index"].append(i_hit)

        # Downstream
        if muFilterHit.GetSystem() == 3 :
            mu["d"][1].append(1.)
        # Upstream
        else :
            mu["d"][1].append(6.)
            
    for i_hit, scifiHit in enumerate(event.Digi_ScifiHits) :
        modules['Scifi'].GetSiPMPosition(scifiHit.GetDetectorID(), a, b)
        scifi["pos"][0].append(a.X())
        scifi["pos"][1].append(a.Y())
        scifi["pos"][2].append(a.Z())

        # 250 mum in x and y directions? maybe?
        scifi["d"][0].append(250e-4)
        scifi["d"][1].append(250e-4)
        # 1.62 mm in z direction?
        scifi["d"][2].append(0.162)
        
        scifi["vert"].append(scifiHit.isVertical())
        scifi["index"].append(i_hit)
        
        scifi["system"].append(0)
    
    for hit_collection in [mu_ds, mu_us, scifi] :
        for key, item in hit_collection.items() :
            hit_collection[key] = np.array(item)

    # Reconstruct muons until there are not enough hits in downstream muon filter
    while True :
        vertical_ds_hits = mu_ds["vert"].sum()
        if vertical_ds_hits < min_hits :
            break
        horizontal_ds_hits = len(mu_ds["vert"]) - vertical_ds_hits
        if horizontal_ds_hits < min_hits :
            break
        # Skip event if there are no scifi hits
        if len(scifi["index"]) == 0 :
            break
        # Skip event if there are no upstream mufilter hits
        if len(mu_us["index"]) == 0 :
            break


        # Hough transform using downstream muon filter only
        ZX = np.dstack([np.tile(mu_ds["pos"][2][mu_ds["vert"]], muon_weight), np.tile(mu_ds["pos"][0][mu_ds["vert"]], muon_weight)])[0]
        d_ZX = np.dstack([np.tile(mu_ds["d"][2][mu_ds["vert"]], muon_weight), np.tile(mu_ds["d"][0][mu_ds["vert"]], muon_weight)])[0]
        ZY = np.dstack([np.tile(mu_ds["pos"][2][~mu_ds["vert"]], muon_weight), np.tile(mu_ds["pos"][1][~mu_ds["vert"]], muon_weight)])[0]
        d_ZY = np.dstack([np.tile(mu_ds["d"][2][~mu_ds["vert"]], muon_weight), np.tile(mu_ds["d"][1][~mu_ds["vert"]], muon_weight)])[0]

        mu_ds_ZX_hough = h_ZX.fit_randomize(ZX, d_ZX, n_random, False)
        mu_ds_ZY_hough = h_ZY.fit_randomize(ZY, d_ZY, n_random, False)
        
        # If fits are unsuccessful, skip to next event
        if not mu_ds_ZX_hough[4] :
            break
        if not mu_ds_ZY_hough[4] :
            break
        
        # And also using SciFi and upstream muon filter
        ZX = np.dstack([np.concatenate([np.tile(mu_ds["pos"][2][mu_ds["vert"]], muon_weight), 
                                        np.tile(mu_us["pos"][2][mu_us["vert"]], muon_weight), 
                                        scifi["pos"][2][scifi["vert"]]]), 
                        np.concatenate([np.tile(mu_ds["pos"][0][mu_ds["vert"]], muon_weight),
                                        np.tile(mu_us["pos"][0][mu_us["vert"]], muon_weight), 
                                        scifi["pos"][0][scifi["vert"]]])])[0]
        d_ZX = np.dstack([np.concatenate([np.tile(mu_ds["d"][2][mu_ds["vert"]], muon_weight), 
                                          np.tile(mu_us["d"][2][mu_us["vert"]], muon_weight), 
                                          scifi["d"][2][scifi["vert"]]]), 
                          np.concatenate([np.tile(mu_ds["d"][0][mu_ds["vert"]], muon_weight), 
                                          np.tile(mu_us["d"][0][mu_us["vert"]], muon_weight), 
                                          scifi["d"][0][scifi["vert"]]])])[0]
        ZY = np.dstack([np.concatenate([np.tile(mu_ds["pos"][2][~mu_ds["vert"]], muon_weight), 
                                        np.tile(mu_us["pos"][2][~mu_us["vert"]], muon_weight), 
                                        scifi["pos"][2][~scifi["vert"]]]), 
                        np.concatenate([np.tile(mu_ds["pos"][1][~mu_ds["vert"]], muon_weight), 
                                        np.tile(mu_us["pos"][1][~mu_us["vert"]], muon_weight), 
                                        scifi["pos"][1][~scifi["vert"]]])])[0]
        d_ZY = np.dstack([np.concatenate([np.tile(mu_ds["d"][2][~mu_ds["vert"]], muon_weight), 
                                          np.tile(mu_us["d"][2][~mu_us["vert"]], muon_weight), 
                                          scifi["d"][2][~scifi["vert"]]]), 
                          np.concatenate([np.tile(mu_ds["d"][1][~mu_ds["vert"]], muon_weight), 
                                          np.tile(mu_us["d"][1][~mu_us["vert"]], muon_weight), 
                                          scifi["d"][1][~scifi["vert"]]])])[0]
        
        all_ZX_hough = h_ZX.fit_randomize(ZX, d_ZX, n_random, False)
        all_ZY_hough = h_ZY.fit_randomize(ZY, d_ZY, n_random, False)

        # Find muon ds hits
        i_ZX = np.tile(mu_ds["index"][mu_ds["vert"]], muon_weight)
        i_ZY = np.tile(mu_ds["index"][~mu_ds["vert"]], muon_weight)
        
        ds_track_hits_ZX_i = np.unique(i_ZX[mu_ds_ZX_hough[3]])
        ds_track_hits_ZY_i = np.unique(i_ZY[mu_ds_ZY_hough[3]])

        # If we're not tracking at least three hits in each plane, stop trying.
        if len(ds_track_hits_ZX_i) < 3 and len(ds_track_hits_ZY_i) < 3 :
            break

        # If track doesn't include at least one scifi hit in each planes, stop trying.
        system_ZX = np.concatenate([np.tile(mu_ds["system"][mu_ds["vert"]], muon_weight),
                                    np.tile(mu_us["system"][mu_us["vert"]], muon_weight),
                                    scifi["system"][scifi["vert"]]])
        system_ZY = np.concatenate([np.tile(mu_ds["system"][~mu_ds["vert"]], muon_weight),
                                    np.tile(mu_us["system"][~mu_us["vert"]], muon_weight),
                                    scifi["system"][~scifi["vert"]]])

        if (system_ZX[all_ZX_hough[3]] == 0).sum() == 0 :
            break
        if (system_ZY[all_ZY_hough[3]] == 0).sum() == 0 :
            break

        this_data.append({"mu_ds" : {"ZX" : mu_ds_ZX_hough, "ZY" : mu_ds_ZY_hough}, 
                          "all" : {"ZX" : all_ZX_hough, "ZY" : all_ZY_hough}})
        

        # Remove muon ds hits and repeat
        # Are we missing some kind of offset for the ZY hits?
        delete_mask = np.isin(mu_ds["index"], np.concatenate([ds_track_hits_ZX_i, ds_track_hits_ZY_i]))

        for key, item in mu_ds.items() :
            if hasattr(item[0], '__len__'):
                mu_ds[key] = item[:,~delete_mask]
            else :
                mu_ds[key] = item[~delete_mask]
#        print("ALL ZX HOUGH ", all_ZX_hough)
#        print("ALL ZY HOUGH ",all_ZY_hough)
#        print("mu_ds", mu_ds)

    print("Muons found: {0}".format(len(this_data)))
    out_data.append(this_data)

with open(args.output_file, "wb") as f_out :
    pickle.dump(out_data, f_out)
