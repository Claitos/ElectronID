#! /usr/bin/python
         
#Import relevant packages
import ROOT, math, array                                
from ROOT import TH1F, TH2F, TMath, TTree, TVector3, TVector2
import uproot as up
import awkward as ak

# Define and open files
infile="Data/pythia8.316-1.0_NC_noRad_ep_10x100_q2_10to100_run000.0000.eicrecon.edm4eic.root" 
ofile=ROOT.TFile.Open("Output_Root/TrackAnalysis_Out.root", "RECREATE")

# Open input file and define branches we want to look at with uproot
events_tree = up.open(infile)["events"]

# Get particle information# Get particle information
mcGenStat = events_tree["MCParticles.generatorStatus"].array()
mcMomX = events_tree["MCParticles.momentum.x"].array() 
mcMomY = events_tree["MCParticles.momentum.y"].array()
mcMomZ = events_tree["MCParticles.momentum.z"].array()
mcCharge = events_tree["MCParticles.charge"].array()
mcPdg = events_tree["MCParticles.PDG"].array()
mcMass = events_tree["MCParticles.mass"].array()

# Get reconstructed track information
tracks = events_tree["ReconstructedParticles"].array()
num_events = len(tracks)

trackMomX = events_tree["ReconstructedParticles.momentum.x"].array()
trackMomY = events_tree["ReconstructedParticles.momentum.y"].array()
trackMomZ = events_tree["ReconstructedParticles.momentum.z"].array()
trackMass = events_tree["ReconstructedParticles.mass"].array()
trackPDG = events_tree["ReconstructedParticles.PDG"].array()

# Get reconstructed cluster information
cluster_links = events_tree["_ReconstructedParticles_clusters.index"].array()

# Get assocations between MCParticles and ReconstructedChargedParticles
recoAssoc = events_tree["_ReconstructedParticleLinks_from.index"].array()
simuAssoc = events_tree["_ReconstructedParticleLinks_to.index"].array()


# Add main analysis loop(s) below
num_tracks = 0
for i in range(len(tracks)):
    num_tracks += len(trackMomX[i])
    for j in range(len(trackMomX[i])):
        #print("Track number: ", j)
        #print("Track momentum (x,y,z): ", trackMomX[i][j], trackMomY[i][j], trackMomZ[i][j])
        # Get the index of the associated MCParticle for this reconstructed track
        mc_index = recoAssoc[i][j]
        if mc_index >= 0:
            pass
            #print("Associated MCParticle index: ", mc_index)
            #print("MCParticle momentum (x,y,z): ", mcMomX[i][mc_index], mcMomY[i][mc_index], mcMomZ[i][mc_index])
            #print("MCParticle PDG code: ", mcPdg[i][mc_index])
        else:
            pass
            #print("No associated MCParticle for this reconstructed track.")

print("Total number of events processed: ", num_events)
print("Total number of reconstructed tracks: ", num_tracks)

# Write output histograms to file below                

# Close files
ofile.Close()    