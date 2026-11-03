import edm4hep                                  # pyright: ignore[reportMissingImports]
from podio import Frame                         # pyright: ignore[reportMissingImports]
from podio.reading import get_reader            # pyright: ignore[reportMissingImports]
import ROOT                                     # pyright: ignore[reportMissingImports]
import XRootD                                   # pyright: ignore[reportMissingImports]
import math               
from ROOT import TH1F, TH2F                     # pyright: ignore[reportMissingImports]
from ROOT.Math import PxPyPzEVector             # pyright: ignore[reportMissingImports]


def file_reader(n_files: int, stream: bool = False, debug: bool = False):

    string_list = []
    
    if stream:
        for i in range(n_files):
            file_path = f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_10to100/pythia8.316-1.0_NC_noRad_ep_10x100_q2_10to100_run000.{i:04d}.eicrecon.edm4eic.root"
            string_list.append(file_path)

    else:
        if n_files < 1:
            raise ValueError("Number of files must be at least 1.")
        elif n_files > 3:
            raise ValueError("Number of files must not exceed 3. I have not downloaded more than 3 files yet:)")

        for i in range(n_files):
            file_path = f"Data/pythia8.316-1.0_NC_noRad_ep_10x100_q2_10to100_run000.{i:04d}.eicrecon.edm4eic.root"
            string_list.append(file_path)

    reader = get_reader(string_list)
    frames = reader.get("events")

    return frames


def create_histograms() -> None:
    # Create histograms
    mass_hist = TH1F("mass_hist", "Mass Distribution; Mass (GeV/c^2); Counts", 100, 0, 5)
    cluster_hist = TH1F("cluster_hist", "Number of Clusters; Number of Clusters; Counts", 5, 0, 5)
    track_hist = TH1F("track_hist", "Number of Tracks; Number of Tracks; Counts", 5, 0, 5)
    pT_hist = TH1F("pT_hist", "Transverse Momentum Distribution; pT (GeV/c); Counts", 100, 0, 10)
    charge_hist = TH1F("charge_hist", "Charge Distribution; Charge; Counts", 5, -2, 2)
    particle_per_event_hist = TH1F("particle_per_event_hist", "Number of Reconstructed Particles per Event; Number of Reconstructed Particles; Counts", 100, 0, 100)

    global histograms
    histograms = {
        "mass_hist": mass_hist,
        "cluster_hist": cluster_hist,
        "track_hist": track_hist,
        "pT_hist": pT_hist,
        "charge_hist": charge_hist,
        "particle_per_event_hist": particle_per_event_hist
    }


def write_histograms(filename) -> None:
    outfile = ROOT.TFile(filename, "RECREATE")

    for name, hist in histograms.items():
        hist.Write(name)

    outfile.Close()
    


def main():

    number_of_files = 2
    frames = file_reader(number_of_files, stream = False, debug = False)

    num_events = len(frames)
    num_reco_particles = 0

    # Create histograms
    create_histograms()

    #  Main analysis loop over events
    for i, frame in enumerate(frames):
        if i % 100 == 0:
            print(f"Processing event {i+1}/{num_events}")

        # Access interesting collections
        mc_particles = frame.get("MCParticles")
        reco_particles = frame.get("ReconstructedParticles")
        reco_assoc = frame.get("ReconstructedParticleAssociations")

        num_reco_particles += len(reco_particles)
        histograms["particle_per_event_hist"].Fill(len(reco_particles))


        # Loop over ReconstructedParticles and fill histograms
        for reco_particle in reco_particles:
            momentum = reco_particle.getMomentum()
            pdg = reco_particle.getPDG()
            charge = reco_particle.getCharge()
            mass = reco_particle.getMass()      
            energy = reco_particle.getEnergy() # excatly the same as : energy = math.sqrt(momentum.x**2 + momentum.y**2 + momentum.z**2 + mass**2)

            mom4 = PxPyPzEVector(momentum.x, momentum.y, momentum.z, energy)
            number_of_clusters = reco_particle.getClusters().size()
            number_of_tracks = reco_particle.getTracks().size()

            histograms["mass_hist"].Fill(mass)
            histograms["cluster_hist"].Fill(number_of_clusters)
            histograms["track_hist"].Fill(number_of_tracks)
            histograms["pT_hist"].Fill(mom4.Pt())
            histograms["charge_hist"].Fill(charge)


            if number_of_clusters > 0 and charge == -1:
                # Access the first cluster associated with the reconstructed particle
                cluster = reco_particle.getClusters()[0]
                cluster_energy = cluster.getEnergy()
                cluster_position = cluster.getPosition()
                #print(f"Cluster: Energy={cluster_energy}, Position=({cluster_position.x}, {cluster_position.y}, {cluster_position.z})")

            else:
                #print("No clusters associated with this reconstructed particle.")
                pass


            #print(f"ReconstructedParticle: PDG={pdg}, Charge={charge}, Momentum=({momentum.x}, {momentum.y}, {momentum.z}), Mass={mass}")




    print("Total number of reconstructed particles: ", num_reco_particles)

    write_histograms("Output_Root/TrackAnalysis_Out_Podio.root")


if __name__ == "__main__":
    main()


