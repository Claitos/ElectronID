import edm4hep                                  # pyright: ignore[reportMissingImports]
from podio import Frame                         # pyright: ignore[reportMissingImports]
from podio.reading import get_reader            # pyright: ignore[reportMissingImports]
import ROOT                                     # pyright: ignore[reportMissingImports]
import XRootD                                   # pyright: ignore[reportMissingImports]
import math               
from ROOT import TH1F, TH2F                     # pyright: ignore[reportMissingImports]
from ROOT.Math import PxPyPzEVector             # pyright: ignore[reportMissingImports]

PDG_CODE_ELECTRON = 11
PDG_CODE_PION = 211
PDG_CODE_KAON = 321
PDG_CODE_PROTON = 2212

def file_reader(n_files: int, stream: bool = False, debug: bool = False) -> list[Frame]:

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
    energy_hist = TH1F("energy_hist", "Energy Distribution; Energy (GeV); Counts", 100, 0, 20)
    energy_cluster_hist = TH1F("energy_cluster_hist", "Cluster Energy Distribution; Cluster Energy (GeV); Counts", 100, 0, 20)
    cluster_hist = TH1F("cluster_hist", "Number of Clusters; Number of Clusters; Counts", 5, 0, 5)
    track_hist = TH1F("track_hist", "Number of Tracks; Number of Tracks; Counts", 5, 0, 5)
    pT_hist = TH1F("pT_hist", "Transverse Momentum Distribution; pT (GeV/c); Counts", 100, 0, 10)
    mom3_hist = TH1F("mom3_hist", "Momentum Distribution; Momentum (GeV/c); Counts", 100, 0, 10)
    charge_hist = TH1F("charge_hist", "Charge Distribution; Charge; Counts", 5, -2, 2)
    particle_per_event_hist = TH1F("particle_per_event_hist", "Number of Reconstructed Particles per Event; Number of Reconstructed Particles; Counts", 100, 0, 100)
    EoverP_hist = TH1F("EoverP_hist", "E/p Distribution; E/p; Counts", 100, 0, 3)
    pid_hist = TH1F("pid_hist", "PID Distribution; PID; Counts", 20, 0, 20)
    pid_goodness_hist = TH1F("pid_goodness_hist", "PID Goodness Distribution; PID Goodness; Counts", 100, 0, 1)
    pdg_hist = TH1F("pdg_hist", "Assigned PDG code (Matching probably?); PDG; Counts", 1000, 0, 2500)
    normalized_electron_likelihood_hist = TH1F("normalized_electron_likelihood_hist", "Normalized Electron Likelihood Distribution; Normalized Electron Likelihood; Counts", 100, -1, 1)

    global histograms
    histograms = {
        "mass_hist": mass_hist,
        "cluster_hist": cluster_hist,
        "track_hist": track_hist,
        "pT_hist": pT_hist,
        "charge_hist": charge_hist,
        "particle_per_event_hist": particle_per_event_hist,
        "energy_hist": energy_hist,
        "EoverP_hist": EoverP_hist,
        "energy_cluster_hist": energy_cluster_hist,
        "pid_hist": pid_hist,
        "pid_goodness_hist": pid_goodness_hist,
        "pdg_hist": pdg_hist,
        "normalized_electron_likelihood_hist": normalized_electron_likelihood_hist,
        "mom3_hist": mom3_hist
    }



def write_histograms(filename: str) -> None:
    outfile = ROOT.TFile(filename, "RECREATE")

    for name, hist in histograms.items():
        hist.Write(name)

    outfile.Close()
    


def main():

    number_of_files = 3
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

            charge = reco_particle.getCharge()
            number_of_clusters = reco_particle.getClusters().size()
            number_of_tracks = reco_particle.getTracks().size()

            histograms["cluster_hist"].Fill(number_of_clusters)
            histograms["track_hist"].Fill(number_of_tracks)
            histograms["charge_hist"].Fill(charge)

            if number_of_clusters > 0 and number_of_tracks > 0 and charge == -1:

                momentum = reco_particle.getMomentum()
                mass = reco_particle.getMass()      
                energy = reco_particle.getEnergy() # excatly the same as : energy = math.sqrt(momentum.x**2 + momentum.y**2 + momentum.z**2 + mass**2)   (I checked it and it is correct so i assume its is energy calculated from the momentum and mass)
                mom4 = PxPyPzEVector(momentum.x, momentum.y, momentum.z, energy)
                
                # Access the first cluster associated with the reconstructed particle
                cluster = reco_particle.getClusters()[0]
                cluster_energy = cluster.getEnergy()
                cluster_position = cluster.getPosition()

                EoverP = cluster_energy / mom4.P()
                
                histograms["mass_hist"].Fill(mass)
                histograms["pT_hist"].Fill(mom4.Pt())
                histograms["mom3_hist"].Fill(mom4.P())
                histograms["energy_hist"].Fill(energy)
                histograms["energy_cluster_hist"].Fill(cluster_energy)
                histograms["EoverP_hist"].Fill(EoverP)

                # Access PID information from the reconstructed particle
                pdg = reco_particle.getPDG()
                histograms["pdg_hist"].Fill(abs(pdg))                # PDG codes assigned by MC matching
                pid_goodness = reco_particle.getGoodnessOfPID()
                histograms["pid_goodness_hist"].Fill(pid_goodness)   # -> is always 0.0 -> no PID used

                pids = reco_particle.getParticleIDs()
                number_of_pid = pids.size()
                histograms["pid_hist"].Fill(number_of_pid)

                max_likelihoods = {"electron": 0.0, "pion": 0.0, "kaon": 0.0, "proton": 0.0}
                for pid in pids:
                    assigned_pdg = abs(pid.getPDG())
                    likelihood = pid.getLikelihood()
                    #print(f"Assigned PDG: {assigned_pdg}, Likelihood: {likelihood}")
                    if assigned_pdg == PDG_CODE_ELECTRON:  # Electron
                        max_likelihoods["electron"] = max(max_likelihoods["electron"], likelihood)
                    elif assigned_pdg == PDG_CODE_PION:  # Pion
                        max_likelihoods["pion"] = max(max_likelihoods["pion"], likelihood)
                    elif assigned_pdg == PDG_CODE_KAON:  # Kaon
                        max_likelihoods["kaon"] = max(max_likelihoods["kaon"], likelihood)
                    elif assigned_pdg == PDG_CODE_PROTON:  # Proton
                        max_likelihoods["proton"] = max(max_likelihoods["proton"], likelihood)

                #print(f"Max Likelihoods for Reconstructed Particle: {max_likelihoods}")
                normalized_electron_likelihood = max_likelihoods["electron"] / sum(max_likelihoods.values()) if sum(max_likelihoods.values()) > 0 else -1.0
                histograms["normalized_electron_likelihood_hist"].Fill(normalized_electron_likelihood)

            else:
                #print("No clusters associated with this reconstructed particle.")
                pass



    print("Total number of reconstructed particles: ", num_reco_particles)
    write_histograms("Output_Root/TrackAnalysis_Out_Podio.root")


if __name__ == "__main__":
    main()


