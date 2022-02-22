import argparse
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import copy

import ROOT
a = ROOT.TVector3()
b = ROOT.TVector3()
# Just a plotting tool

def getTruth(event) :
    mc_truth = {}
    
#    for point_collections in [event.MuFilterPoint, event.vetoPoint, event.ScifiPoint ] :
    for point_collections in [event.MuFilterPoint, event.ScifiPoint ] :
        for point in point_collections :
            if point.GetTrackID() < 0 :
                continue
            if point.GetTrackID() not in mc_truth :
                mc_truth[point.GetTrackID()] = [point.PdgCode(), [], [], [], [], [], [], []]
            else :
                mc_truth[point.GetTrackID()][1].append(point.GetX())
                mc_truth[point.GetTrackID()][2].append(point.GetY())
                mc_truth[point.GetTrackID()][3].append(point.GetZ())
                mc_truth[point.GetTrackID()][4].append(point.GetTime())
                mc_truth[point.GetTrackID()][5].append(point.GetPx())
                mc_truth[point.GetTrackID()][6].append(point.GetPy())
                mc_truth[point.GetTrackID()][7].append(point.GetPz())
    return mc_truth

mc_colors = {}
mc_colors[13] = "tab:red"
mc_colors[22] = "tab:olive"
mc_colors[11] = "tab:olive"

def draw_truth(mc_truth, dims) :
    for trackID, track in mc_truth.items() :
        if abs(track[0]) in mc_colors :
            this_color = mc_colors[abs(track[0])]
        else :
            this_color = "tab:purple"
        plt.plot(track[1+dims[0]], track[1+dims[1]], color = this_color)

def draw_hits(ax, x, d, color, label) :
    for i in range(len(x)) :
        ax.add_patch(Rectangle(
            xy=(x[i][0] - d[i][0]/2, x[i][1] - d[i][1]/2.) ,width=d[i][0], height=d[i][1],
            linewidth=0, color=color, fill=True))
        if not i and label is not None :
            handles, labels = ax.get_legend_handles_labels()
            patch = mpatches.Patch(color=color, label=label)
            handles.append(patch) 

# Geometry stuff
from rootpyPickler import Unpickler

parser = argparse.ArgumentParser()
parser.add_argument("input_digi_file")
parser.add_argument("geo_file")
parser.add_argument("event_number", type = int)
parser.add_argument("--reco_file")
parser.add_argument("--out_name")
parser.add_argument("--truth", type = bool, default = False)
parser.add_argument("--h6", dest = "h6", action = "store_true")
parser.set_defaults(h6 = False)

args = parser.parse_args()

fgeo = ROOT.TFile.Open(args.geo_file)
upkl    = Unpickler(fgeo)
snd_geo = upkl.load('ShipGeo')
run = ROOT.FairRunSim()

import shipLHC_conf as sndDet_conf
modules = sndDet_conf.configure(run,snd_geo)
sGeo = fgeo.FAIRGeom
modules['Scifi'].SiPMmapping()

treeName = "cbmsim"
if args.h6 :
    treeName = "rawConv"

event = ROOT.TChain(treeName)

event.Add(args.input_digi_file)

print(event.GetEntries())
print(args.event_number)

if args.reco_file is not None :
    reco_event = ROOT.TChain(treeName)
    reco_event.Add(args.reco_file)
    event.AddFriend(reco_event)

event.GetEntry(args.event_number)

# Read hits
# For downstream muon filter hits
base_hit_dict = {"pos" : [[], [], []], 
                 "d" : [[], [], []], 
                 "vert" : [], 
                 "index" : [],
                 "system" : [], 
                 "in_track" : []}

mu_ds = copy.deepcopy(base_hit_dict)

# For upstream muon filter hits
mu_us = copy.deepcopy(base_hit_dict)

# For scifi hits
scifi = copy.deepcopy(base_hit_dict)

# For veto hits
veto = copy.deepcopy(base_hit_dict)

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
    mu["pos"][0].append((a+b).X()/2.)
    mu["pos"][1].append((a+b).Y()/2.)
    mu["pos"][2].append((a+b).Z()/2.)
    
    mu["vert"].append(muFilterHit.isVertical())
    mu["system"].append(muFilterHit.GetSystem())
        
    mu["d"][0].append(1.+np.abs((a-b).X()))
    mu["d"][2].append(1.+np.abs((a-b).Z()))
        
    mu["index"].append(i_hit)
    
    # Downstream
    if muFilterHit.GetSystem() == 3 :
        mu["d"][1].append(1.+np.abs((a-b).Y()))
    # Upstream or veto
    else :
        mu["d"][1].append(6.+np.abs((a-b).Y()))

    in_track = False
#    if args.reco_file is not None :
#        for reco_track in event.Reco_MuonTracks :
#            if muFilterHit.GetDetectorID() in reco_track.getHits() :
#                in_track = True

    mu["in_track"].append(in_track)
            
for i_hit, scifiHit in enumerate(event.Digi_ScifiHits) :
    modules['Scifi'].GetSiPMPosition(scifiHit.GetDetectorID(), a, b)
    scifi["pos"][0].append((a+b).X()/2.)
    scifi["pos"][1].append((a+b).Y()/2.)
    scifi["pos"][2].append((a+b).Z()/2.)

    # 250 mum in x and y directions? maybe?
    scifi["d"][0].append(250e-4+np.abs((a-b).X()))
    scifi["d"][1].append(250e-4+np.abs((a-b).Y()))
    # 1.62 mm in z direction?
    scifi["d"][2].append(0.162+np.abs((a-b).Z()))
    
    scifi["vert"].append(scifiHit.isVertical())
    scifi["index"].append(i_hit)
        
    scifi["system"].append(0)

    in_track = False
#    if args.reco_file is not None :
#        for reco_track in event.Reco_MuonTracks :
#            if scifiHit.GetDetectorID() in reco_track.getHits() :
#                in_track = True

    scifi["in_track"].append(in_track)

print(scifi["pos"][0])
print(scifi["pos"][1])
print(scifi["pos"][2])
    
for hit_collection in [mu_ds, mu_us, veto, scifi] :
    for key, item in hit_collection.items() :
        hit_collection[key] = np.array(item)

plt.figure(figsize = [1.5*6.4, 1.5*4.8])
ax_zx = plt.subplot(2, 1, 1)
#if args.truth is not None :
#    if args.truth :
#        for trackID, track in truth.items() :
#            plt.plot([p[2] for p in track[1]], [p[0] for p in track[1]], linewidth = 1., color = truth_ana.track_color(track[0]), alpha = 0.5)

if len(mu_ds["vert"]) :
    draw_hits(ax_zx, 
              np.dstack([mu_ds["pos"][2][mu_ds["vert"]], mu_ds["pos"][0][mu_ds["vert"]]])[0], 
              np.dstack([mu_ds["d"][2][mu_ds["vert"]], mu_ds["d"][0][mu_ds["vert"]]])[0],
              "tab:blue", "MuFilter")
if len(scifi["vert"]) :
    draw_hits(ax_zx, 
              np.dstack([scifi["pos"][2][scifi["vert"]], scifi["pos"][0][scifi["vert"]]])[0], 
              np.dstack([scifi["d"][2][scifi["vert"]], scifi["d"][0][scifi["vert"]]])[0],
              "tab:red", "Scifi")

if args.reco_file is not None :
    if len(mu_ds["vert"]) :
        draw_hits(ax_zx,
                  np.dstack([mu_ds["pos"][2][np.logical_and(mu_ds["vert"], mu_ds["in_track"])], mu_ds["pos"][0][np.logical_and(mu_ds["vert"], mu_ds["in_track"])]])[0], 
                  np.dstack([mu_ds["d"][2][np.logical_and(mu_ds["vert"], mu_ds["in_track"])], mu_ds["d"][0][np.logical_and(mu_ds["vert"], mu_ds["in_track"])]])[0],
                  "tab:orange", "MuFilter")

    if len(scifi["vert"]) :
        draw_hits(ax_zx, 
                  np.dstack([scifi["pos"][2][np.logical_and(scifi["vert"], scifi["in_track"])], scifi["pos"][0][np.logical_and(scifi["vert"], scifi["in_track"])]])[0], 
                  np.dstack([scifi["d"][2][np.logical_and(scifi["vert"], scifi["in_track"])], scifi["d"][0][np.logical_and(scifi["vert"], scifi["in_track"])]])[0],
                  "tab:orange", "Scifi")
                  
if args.truth :
    mc_truth = getTruth(event)
    draw_truth(mc_truth, [2, 0])

#plt.xlabel("z [cm]")
plt.ylabel("x [cm]")


ax_zy = plt.subplot(2, 1, 2)
#if args.truth is not None :
#    if args.truth :
#        for trackID, track in truth.items() :
#            plt.plot([p[2] for p in track[1]], [p[1] for p in track[1]], linewidth = 1., color = truth_ana.track_color(track[0]), alpha = 0.5)

if len(mu_ds["vert"]) :
    draw_hits(ax_zy, 
              np.dstack([mu_ds["pos"][2][~mu_ds["vert"]], mu_ds["pos"][1][~mu_ds["vert"]]])[0], 
              np.dstack([mu_ds["d"][2][~mu_ds["vert"]], mu_ds["d"][1][~mu_ds["vert"]]])[0],
              "tab:blue", "MuFilter")
if len(veto["vert"]) :
    draw_hits(ax_zy, 
              np.dstack([veto["pos"][2][~veto["vert"]], veto["pos"][1][~veto["vert"]]])[0], 
              np.dstack([veto["d"][2][~veto["vert"]], veto["d"][1][~veto["vert"]]])[0],
              "tab:brown", None)
if len(mu_us["vert"]) :
    draw_hits(ax_zy, 
              np.dstack([mu_us["pos"][2][~mu_us["vert"]], mu_us["pos"][1][~mu_us["vert"]]])[0], 
              np.dstack([mu_us["d"][2][~mu_us["vert"]], mu_us["d"][1][~mu_us["vert"]]])[0],
              "tab:blue", None)
if len(scifi["vert"]) :
    draw_hits(ax_zy, 
              np.dstack([scifi["pos"][2][~scifi["vert"]], scifi["pos"][1][~scifi["vert"]]])[0], 
              np.dstack([scifi["d"][2][~scifi["vert"]], scifi["d"][1][~scifi["vert"]]])[0],
              "tab:red", "Scifi")

if args.reco_file is not None :
    if len(mu_ds["vert"]) :
        draw_hits(ax_zy, 
                  np.dstack([mu_ds["pos"][2][np.logical_and(~mu_ds["vert"], mu_ds["in_track"])], mu_ds["pos"][1][np.logical_and(~mu_ds["vert"], mu_ds["in_track"])]])[0], 
                  np.dstack([mu_ds["d"][2][np.logical_and(~mu_ds["vert"], mu_ds["in_track"])], mu_ds["d"][1][np.logical_and(~mu_ds["vert"], mu_ds["in_track"])]])[0],
                  "tab:orange", "MuFilter")
    if len(veto["vert"]) :
        draw_hits(ax_zy, 
                  np.dstack([veto["pos"][2][np.logical_and(~veto["vert"], veto["in_track"])], veto["pos"][1][np.logical_and(~veto["vert"], veto["in_track"])]])[0], 
                  np.dstack([veto["d"][2][np.logical_and(~veto["vert"], veto["in_track"])], veto["d"][1][np.logical_and(~veto["vert"], veto["in_track"])]])[0],
                  "tab:orange", None)
    if len(mu_us["vert"]) :
        draw_hits(ax_zy, 
                  np.dstack([mu_us["pos"][2][np.logical_and(~mu_us["vert"], mu_us["in_track"])], mu_us["pos"][1][np.logical_and(~mu_us["vert"], mu_us["in_track"])]])[0], 
                  np.dstack([mu_us["d"][2][np.logical_and(~mu_us["vert"], mu_us["in_track"])], mu_us["d"][1][np.logical_and(~mu_us["vert"], mu_us["in_track"])]])[0],
                  "tab:orange", None)
    if len(scifi["vert"]) :
        draw_hits(ax_zy, 
                  np.dstack([scifi["pos"][2][np.logical_and(~scifi["vert"], scifi["in_track"])], scifi["pos"][1][np.logical_and(~scifi["vert"], scifi["in_track"])]])[0], 
                  np.dstack([scifi["d"][2][np.logical_and(~scifi["vert"], scifi["in_track"])], scifi["d"][1][np.logical_and(~scifi["vert"], scifi["in_track"])]])[0],
                  "tab:orange", "Scifi")

if args.truth :
    draw_truth(mc_truth, [2, 1])


if args.reco_file is not None :
    for reco_track in event.Reco_KalmanTracks :
        points_z = []
        points_x = []
        points_y = []
        for i in range(reco_track.getNumPointsWithMeasurement()):
            state = reco_track.getFittedState(i)
            pos    = state.getPos()
            points_z.append(pos[2])
            points_x.append(pos[0])
            points_y.append(pos[1])
        ax_zx.plot(points_z, points_x, color = "tab:green")
        ax_zy.plot(points_z, points_y, color = "tab:green")


#    for reco_track in event.Reco_MuonTracks :
#        ax_zx.plot([reco_track.getStart().Z(), reco_track.getStop().Z()], [reco_track.getStart().X(), reco_track.getStop().X()], color = "tab:red")
#        ax_zy.plot([reco_track.getStart().Z(), reco_track.getStop().Z()], [reco_track.getStart().Y(), reco_track.getStop().Y()], color = "tab:red")

#ax_zx.legend()

plt.xlabel("z [cm]")
plt.ylabel("y [cm]")

if args.h6 :
    z_lims = [-60, 380]
else :
    z_lims = [340+-50, 340+250]


ax_zx.set_xlim(z_lims[0], z_lims[1])
ax_zx.set_ylim(-80, 0)

ax_zy.set_xlim(z_lims[0], z_lims[1])
ax_zy.set_ylim(0, 80)

plt.tight_layout()

if args.out_name is not None :
    plt.savefig(args.out_name)

plt.show()
