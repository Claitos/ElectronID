import edm4hep                                  # pyright: ignore[reportMissingImports]
from edm4hep import ReconstructedParticle       # pyright: ignore[reportMissingImports]
import math
import ROOT                                     # pyright: ignore[reportMissingImports]             
from ROOT import TH1F, TH2F                     # pyright: ignore[reportMissingImports]
from ROOT.Math import PxPyPzEVector             # pyright: ignore[reportMissingImports]
ROOT.TH1.AddDirectory(False)
from Libs.particle_calculations import calculate_elctron_likelihood, calculate_isolation_energy
from Libs.file_reader import file_reader


def create_histograms() -> None:

    global signal_histograms, background_histograms
    nbins = 500
    pT_range = (0, 20)
    EoverP_range = (0, 3)

    signal_histograms = {
        "pT_hist": TH1F("pT_hist", "Transverse Momentum Distribution; pT (GeV/c); Counts", nbins, *pT_range),
        "EoverP_hist": TH1F("EoverP_hist", "E/p Distribution; E/p; Counts", nbins, *EoverP_range),
        "normalized_electron_likelihood_hist": TH1F("normalized_electron_likelihood_hist", "Normalized Electron Likelihood Distribution; Normalized Electron Likelihood; Counts", nbins, 0, 1.1),
        "isolation_energy_ratio_hist": TH1F("isolation_energy_ratio_hist", "Isolation Energy Ratio Distribution; Isolation Energy Ratio; Counts", nbins, 0, 1.1),
    }

    background_histograms = {
        "pT_hist": TH1F("pT_hist", "Transverse Momentum Distribution; pT (GeV/c); Counts", nbins, *pT_range),
        "EoverP_hist": TH1F("EoverP_hist", "E/p Distribution; E/p; Counts", nbins, *EoverP_range),
        "normalized_electron_likelihood_hist": TH1F("normalized_electron_likelihood_hist", "Normalized Electron Likelihood Distribution; Normalized Electron Likelihood; Counts", nbins, 0, 1.1),
        "isolation_energy_ratio_hist": TH1F("isolation_energy_ratio_hist", "Isolation Energy Ratio Distribution; Isolation Energy Ratio; Counts", nbins, 0, 1.1),
    }

    global signal_pdfs, background_pdfs
    signal_pdfs = {
        "pT_hist": TH1F("pT_hist", "Transverse Momentum Distribution; pT (GeV/c); Probability Density", nbins, *pT_range),
        "EoverP_hist": TH1F("EoverP_hist", "E/p Distribution; E/p; Probability Density", nbins, *EoverP_range),
        "normalized_electron_likelihood_hist": TH1F("normalized_electron_likelihood_hist", "Normalized Electron Likelihood Distribution; Normalized Electron Likelihood; Probability Density", nbins, 0, 1.1),
        "isolation_energy_ratio_hist": TH1F("isolation_energy_ratio_hist", "Isolation Energy Ratio Distribution; Isolation Energy Ratio; Probability Density", nbins, 0, 1.1),
    }
    background_pdfs = {
        "pT_hist": TH1F("pT_hist", "Transverse Momentum Distribution; pT (GeV/c); Probability Density", nbins, *pT_range),
        "EoverP_hist": TH1F("EoverP_hist", "E/p Distribution; E/p; Probability Density", nbins, *EoverP_range),
        "normalized_electron_likelihood_hist": TH1F("normalized_electron_likelihood_hist", "Normalized Electron Likelihood Distribution; Normalized Electron Likelihood; Probability Density", nbins, 0, 1.1),
        "isolation_energy_ratio_hist": TH1F("isolation_energy_ratio_hist", "Isolation Energy Ratio Distribution; Isolation Energy Ratio; Probability Density", nbins, 0, 1.1),
    }

    
def normalize_histograms() -> None:

    for hist in signal_pdfs.values():

        integral = hist.Integral()

        if integral > 0:
            hist.Scale(1.0 / integral, "width")

    for hist in background_pdfs.values():

        integral = hist.Integral()

        if integral > 0:
            hist.Scale(1.0 / integral, "width")


def check_normalization() -> None:

    print("\nChecking histogram normalization:")

    for category, histograms in [
        ("Signal", signal_pdfs),
        ("Background", background_pdfs)
    ]:

        for name, hist in histograms.items():

            integral = hist.Integral("width")

            print(f"{category:10s} | {name:40s} | Integral = {integral:.6f}")


def write_histograms(filename: str) -> None:
    outfile = ROOT.TFile(filename, "RECREATE")

    # Create directories
    signal_dir = outfile.mkdir("Signal")
    background_dir = outfile.mkdir("Background")

    signal_counts_dir = signal_dir.mkdir("Counts")
    signal_pdfs_dir = signal_dir.mkdir("PDFs")

    background_counts_dir = background_dir.mkdir("Counts")
    background_pdfs_dir = background_dir.mkdir("PDFs")

    # Write signal histograms
    signal_counts_dir.cd()
    for name, hist in signal_histograms.items():
        hist.Write(name)

    signal_pdfs_dir.cd()
    for name, hist in signal_pdfs.items():
        hist.Write(name)

    # Write background histograms
    background_counts_dir.cd()
    for name, hist in background_histograms.items():
        hist.Write(name)

    background_pdfs_dir.cd()
    for name, hist in background_pdfs.items():
        hist.Write(name)

    outfile.Close()





def main():

    number_of_files = 100
    frames = file_reader(number_of_files, stream = True, debug = False)

    num_events = len(frames)
    num_reco_particles = 0
    num_scattered_electrons_total = 0
    num_background_particles_total = 0

    # Create histograms
    create_histograms()

    #  Main analysis loop over events
    for i, frame in enumerate(frames):
        if i % 100 == 0:
            print(f"Processing event {i+1}/{num_events}")

        # Access interesting collections
        # mc_particles = frame.get("MCParticles")
        # reco_assoc = frame.get("ReconstructedParticleAssociations")

        reco_particles = frame.get("ReconstructedParticles")
        scattered_electrons_Truth = frame.get("ScatteredElectronsTruth")

        num_reco_particles += len(reco_particles)
        num_scattered_electrons = len(scattered_electrons_Truth)

        if num_scattered_electrons == 0:   # If there are no scattered electrons in the event, skip to the next event
            continue

        # Loop over ReconstructedParticles and fill histograms
        for reco_particle in reco_particles:

            charge = reco_particle.getCharge()
            number_of_clusters = reco_particle.getClusters().size()
            number_of_tracks = reco_particle.getTracks().size()

            if number_of_clusters == 0 or number_of_tracks == 0:
                continue

            if charge != -1:
                continue

            momentum = reco_particle.getMomentum()   
            energy = reco_particle.getEnergy() # excatly the same as : energy = math.sqrt(momentum.x**2 + momentum.y**2 + momentum.z**2 + mass**2)   (I checked it and it is correct so i assume its is energy calculated from the momentum and mass)
            mom4 = PxPyPzEVector(momentum.x, momentum.y, momentum.z, energy) 
            
            # Access the first cluster associated with the reconstructed particle  (Check in cluster_hist that there is either 0 or 1 cluster associated with the reconstructed particle)
            cluster = reco_particle.getClusters()[0]
            cluster_energy = cluster.getEnergy()
            cluster_position = cluster.getPosition()

            EoverP = cluster_energy / mom4.P()
            isolation_energy = calculate_isolation_energy(cluster_position, reco_particles)


            # Access PID information from the reconstructed particle
            pids = reco_particle.getParticleIDs()
            number_of_pid = pids.size()

            pid_flag = True   #   Before checking if there is PID information, we assume that there is PID information.  If there is no PID information, we set the flag to False and skip the likelihood calculation.
            if number_of_pid > 0:
                normalized_electron_likelihood = calculate_elctron_likelihood(reco_particle)
            else:
                # Here i have to careful since no PID carries different meaning then likelyhood == 0.
                pid_flag = False


            # Fill histograms based on being scattered electron or not
            esct_id = scattered_electrons_Truth[0].getObjectID().index
            reco_id = reco_particle.getObjectID().index

            if esct_id == reco_id:
                #print(f"Scattered electron is matched to a reconstructed particle with PDG {reco_particle.getPDG()}.  {scattered_electrons_Truth[0].getPDG() == reco_particle.getPDG()}")
                num_scattered_electrons_total += 1
                signal_histograms["pT_hist"].Fill(mom4.Pt())
                signal_histograms["EoverP_hist"].Fill(EoverP)
                signal_histograms["isolation_energy_ratio_hist"].Fill(cluster_energy / isolation_energy)

                signal_pdfs["pT_hist"].Fill(mom4.Pt())
                signal_pdfs["EoverP_hist"].Fill(EoverP)
                signal_pdfs["isolation_energy_ratio_hist"].Fill(cluster_energy / isolation_energy)
                
                if pid_flag:
                    signal_histograms["normalized_electron_likelihood_hist"].Fill(normalized_electron_likelihood)
                    signal_pdfs["normalized_electron_likelihood_hist"].Fill(normalized_electron_likelihood)

            else:
                num_background_particles_total += 1
                background_histograms["pT_hist"].Fill(mom4.Pt())
                background_histograms["EoverP_hist"].Fill(EoverP)
                background_histograms["isolation_energy_ratio_hist"].Fill(cluster_energy / isolation_energy)

                background_pdfs["pT_hist"].Fill(mom4.Pt())
                background_pdfs["EoverP_hist"].Fill(EoverP)
                background_pdfs["isolation_energy_ratio_hist"].Fill(cluster_energy / isolation_energy)

                if pid_flag:
                    background_histograms["normalized_electron_likelihood_hist"].Fill(normalized_electron_likelihood)
                    background_pdfs["normalized_electron_likelihood_hist"].Fill(normalized_electron_likelihood)



    print("\nTotal number of reconstructed particles: ", num_reco_particles)
    print("Total number of scattered electrons matched to reconstructed particles: ", num_scattered_electrons_total)
    print("Total number of background particles in events with scattered electron: ", num_background_particles_total)
    normalize_histograms()
    check_normalization()
    write_histograms("Output_Root/eID_MC_hists.root")


if __name__ == "__main__":
    main()


