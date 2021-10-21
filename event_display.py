import argparse
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches

import ROOT
a = ROOT.TVector3()
b = ROOT.TVector3()
# Just a plotting tool

def draw_hits(ax, x, d, color, label) :
    for i in range(len(x)) :
        ax.add_patch(Rectangle(
            xy=(x[i][0] - d[i][0]/2, x[i][1] - d[i][1]/2.) ,width=d[i][0], height=d[i][1],
            linewidth=0, color=color, fill=True))
        if not i and label is not None :
            handles, labels = ax.get_legend_handles_labels()
            patch = mpatches.Patch(color=color, label=label)
            handles.append(patch) 

#from hough_transform import hough
import truth_ana

# Geometry stuff
from rootpyPickler import Unpickler

parser = argparse.ArgumentParser()
parser.add_argument("input_digi_file")
parser.add_argument("geo_file")
parser.add_argument("event_number", type = int)
parser.add_argument("--reco_file")
parser.add_argument("--out_name")
parser.add_argument("--truth", type = bool)

args = parser.parse_args()

fgeo = ROOT.TFile.Open(args.geo_file)
upkl    = Unpickler(fgeo)
snd_geo = upkl.load('ShipGeo')
run = ROOT.FairRunSim()

import shipLHC_conf as sndDet_conf
modules = sndDet_conf.configure(run,snd_geo)
sGeo = fgeo.FAIRGeom
modules['Scifi'].SiPMmapping()

event = ROOT.TChain("cbmsim")

event.Add(args.input_digi_file)

event.GetEntry(args.event_number)

truth = truth_ana.getTrueTracks(event)

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

# For veto hits
veto = {"pos" : [[], [], []], 
        "d" : [[], [], []], 
        "vert" : [],
        "index" : [], 
        "system" : []}

# Loop through hits
for i_hit, muFilterHit in enumerate(event.Digi_MuFilterHits) :
    
    # Don't use veto for fitting
    if muFilterHit.GetSystem() == 1 :
        mu = veto
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
    # Upstream or veto
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
    
for hit_collection in [mu_ds, mu_us, veto, scifi] :
    for key, item in hit_collection.items() :
        hit_collection[key] = np.array(item)

plt.figure(figsize = [1.5*6.4, 1.5*4.8])
ax_zx = plt.subplot(2, 1, 1)
if args.truth is not None :
    if args.truth :
        for trackID, track in truth.items() :
            plt.plot([p[2] for p in track[1]], [p[0] for p in track[1]], linewidth = 1., color = truth_ana.track_color(track[0]), alpha = 0.5)

draw_hits(ax_zx, 
          np.dstack([mu_ds["pos"][2][mu_ds["vert"]], mu_ds["pos"][0][mu_ds["vert"]]])[0], 
          np.dstack([mu_ds["d"][2][mu_ds["vert"]], mu_ds["d"][0][mu_ds["vert"]]])[0],
          "tab:blue", "MuFilter")
draw_hits(ax_zx, 
          np.dstack([scifi["pos"][2][scifi["vert"]], scifi["pos"][0][scifi["vert"]]])[0], 
          np.dstack([scifi["d"][2][scifi["vert"]], scifi["d"][0][scifi["vert"]]])[0],
          "tab:cyan", "Scifi")

#plt.xlabel("z [cm]")
plt.ylabel("x [cm]")


ax_zy = plt.subplot(2, 1, 2)
if args.truth is not None :
    if args.truth :
        for trackID, track in truth.items() :
            plt.plot([p[2] for p in track[1]], [p[1] for p in track[1]], linewidth = 1., color = truth_ana.track_color(track[0]), alpha = 0.5)

draw_hits(ax_zy, 
          np.dstack([mu_ds["pos"][2][~mu_ds["vert"]], mu_ds["pos"][1][~mu_ds["vert"]]])[0], 
          np.dstack([mu_ds["d"][2][~mu_ds["vert"]], mu_ds["d"][1][~mu_ds["vert"]]])[0],
          "tab:blue", "MuFilter")
draw_hits(ax_zy, 
          np.dstack([mu_us["pos"][2][~mu_us["vert"]], mu_us["pos"][1][~mu_us["vert"]]])[0], 
          np.dstack([mu_us["d"][2][~mu_us["vert"]], mu_us["d"][1][~mu_us["vert"]]])[0],
          "tab:blue", None)
draw_hits(ax_zy, 
          np.dstack([veto["pos"][2][~veto["vert"]], veto["pos"][1][~veto["vert"]]])[0], 
          np.dstack([veto["d"][2][~veto["vert"]], veto["d"][1][~veto["vert"]]])[0],
          "tab:brown", None)
draw_hits(ax_zy, 
          np.dstack([scifi["pos"][2][~scifi["vert"]], scifi["pos"][1][~scifi["vert"]]])[0], 
          np.dstack([scifi["d"][2][~scifi["vert"]], scifi["d"][1][~scifi["vert"]]])[0],
          "tab:cyan", "Scifi")

if args.reco_file is not None :
    with open(args.reco_file, "rb") as f :
        reco = pickle.load(f)
    for i_track, track in enumerate(reco[args.event_number]) :
        if not i_track :
            label = "Hough track"
        else :
            label = None
        ax_zx.plot(track['all']['ZX'][2][:,0], track['all']['ZX'][2][:,1], label = label, color = "tab:green", alpha = 0.5) 
        ax_zy.plot(track['all']['ZY'][2][:,0], track['all']['ZY'][2][:,1], color = "tab:green", alpha = 0.5) 
        print(track['all']['ZX'][2])
        print(track['all']['ZY'][2])

ax_zx.legend()

plt.xlabel("z [cm]")
plt.ylabel("y [cm]")

ax_zx.set_xlim(-50, 250)
ax_zx.set_ylim(-80, 0)

ax_zy.set_xlim(-50, 250)
ax_zy.set_ylim(0, 80)
plt.tight_layout()

if args.out_name is not None :
    plt.savefig(args.out_name)

plt.show()
