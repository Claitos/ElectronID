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
    # Create histograms
    pT_hist = TH1F("pT_hist", "Transverse Momentum Distribution; pT (GeV/c); Counts", 100, 0, 10)
    EoverP_hist = TH1F("EoverP_hist", "E/p Distribution; E/p; Counts", 100, 0, 3)
    normalized_electron_likelihood_hist = TH1F("normalized_electron_likelihood_hist", "Normalized Electron Likelihood Distribution; Normalized Electron Likelihood; Counts", 100, 0, 1.1)
    isolation_energy_ratio_hist = TH1F("isolation_energy_ratio_hist", "Isolation Energy Ratio Distribution; Isolation Energy Ratio; Counts", 100, 0, 1.1)

    p_signal_hist = TH1F("p_signal_hist", "Signal Probability Distribution; Signal Probability; Counts", 100, 0, 1.1)

    global histograms
    histograms = {
        "pT_hist": pT_hist,
        "EoverP_hist": EoverP_hist,
        "normalized_electron_likelihood_hist": normalized_electron_likelihood_hist,
        "isolation_energy_ratio_hist": isolation_energy_ratio_hist,
        "p_signal_hist": p_signal_hist
    }


def write_histograms(filename: str) -> None:
    outfile = ROOT.TFile(filename, "RECREATE")

    for name, hist in histograms.items():
        hist.Write(name)

    outfile.Close()





def evaluate_pdf(pdf: TH1F, value: float) -> float:

    bin_number = pdf.FindBin(value)

    density = pdf.GetBinContent(bin_number)

    # Small relative density floor
    density_floor = max(pdf.GetMaximum() * 1e-12, 1e-300)

    return max(density, density_floor)


def calculate_log_likelihoods(observables: dict) -> tuple[float, float]:

    log_likelihood_signal = 0.0
    log_likelihood_background = 0.0

    for name, value in observables.items():

        signal_pdf = signal_pdfs[name]
        background_pdf = background_pdfs[name]

        # Get probability densities
        p_signal = evaluate_pdf(signal_pdf, value)
        p_background = evaluate_pdf(background_pdf, value)

        log_likelihood_signal += math.log(p_signal)
        log_likelihood_background += math.log(p_background)

    return log_likelihood_signal, log_likelihood_background



def calculate_signal_probability(log_likelihood_signal: float, log_likelihood_background: float, prior_signal: float, prior_background: float) -> float:

    log_odds = (log_likelihood_signal - log_likelihood_background + math.log(prior_signal / prior_background))

    return 1.0 / (1.0 + math.exp(-log_odds))



def calculate_priors() -> tuple[float, float]:
    
    # Calculate priors based on the number of signal and background events in the MC sample
    num_signal_events = signal_pdfs["pT_hist"].GetEntries()
    num_background_events = background_pdfs["pT_hist"].GetEntries()

    total_events = num_signal_events + num_background_events

    prior_signal = num_signal_events / total_events
    prior_background = num_background_events / total_events

    return prior_signal, prior_background



def main():

    number_of_files = 3
    frames = file_reader(number_of_files, stream = False, debug = False)

    num_events = len(frames)
    num_reco_particles = 0
    num_escattered_electron_candidates = 0   # above 0.9 signal probability

    # Create histograms
    create_histograms()

    # Read in MC pdfs
    infile = ROOT.TFile("Output_Root/eID_MC_hists.root", "READ")
    global signal_pdfs, background_pdfs
    signal_pdfs = infile.Get("Signal/PDFs")
    background_pdfs = infile.Get("Background/PDFs")

    # Calculate priors based on the number of signal and background events in the MC sample
    prior_signal, prior_background = calculate_priors()
    print(f"Prior signal: {prior_signal}, Prior background: {prior_background} \n")

    #  Main analysis loop over events
    for i, frame in enumerate(frames):
        if i % 100 == 0:
            print(f"Processing event {i+1}/{num_events}")

        # Access interesting collections
        mc_particles = frame.get("MCParticles")
        reco_particles = frame.get("ReconstructedParticles")
        reco_assoc = frame.get("ReconstructedParticleAssociations")

        num_reco_particles += len(reco_particles)

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

            pid_flag = True     # before checking if there is PID information, we assume that there is PID information.
            if number_of_pid > 0:
                normalized_electron_likelihood = calculate_elctron_likelihood(reco_particle)
            else:
                # Here i have to careful since no PID carries different meaning then likelyhood == 0.
                pid_flag = False

            histograms["pT_hist"].Fill(mom4.Pt())
            histograms["EoverP_hist"].Fill(EoverP)
            histograms["isolation_energy_ratio_hist"].Fill(cluster_energy / isolation_energy)
            histograms["normalized_electron_likelihood_hist"].Fill(normalized_electron_likelihood)

            if pid_flag:
                observables = {"pT_hist": mom4.Pt(),
                                "EoverP_hist": EoverP,
                                "normalized_electron_likelihood_hist": normalized_electron_likelihood,
                                "isolation_energy_ratio_hist": cluster_energy / isolation_energy
                                }
            else:
                observables = {"pT_hist": mom4.Pt(),
                                "EoverP_hist": EoverP,
                                "isolation_energy_ratio_hist": cluster_energy / isolation_energy
                                }

            log_sig, log_bkg = calculate_log_likelihoods(observables)
            p_signal = calculate_signal_probability(log_sig, log_bkg, prior_signal, prior_background)

            histograms["p_signal_hist"].Fill(p_signal)

            if p_signal > 0.9:
                num_escattered_electron_candidates += 1



    print("\nTotal number of reconstructed particles: ", num_reco_particles)
    print("Total number of scattered electron candidates (p_signal > 0.9): ", num_escattered_electron_candidates)
    write_histograms("Output_Root/eID_Rec_hists.root")


if __name__ == "__main__":
    main()


