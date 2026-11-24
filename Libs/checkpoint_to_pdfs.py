from pathlib import Path
from ROOT import TH1F                   # pyright: ignore[reportMissingImports]
import ROOT                           # pyright: ignore[reportMissingImports]


## ---------------------------------------------------------------------------
OUTPUT_DIR = Path("Output_Root_Checkpointed")
FINAL_OUTPUT_FILE = OUTPUT_DIR / "eID_MC_hists.root"
## ---------------------------------------------------------------------------

def read_histogram_set(root_file, directory: str, names: list[str]) -> dict[str, TH1F]:
    histograms: dict[str, TH1F] = {}

    for name in names:
        stored_hist = root_file.Get(f"{directory}/Counts/{name}")

        # Clone so the histograms remain valid after closing the ROOT input file.
        hist = stored_hist.Clone(f"{directory.lower()}_{name}")
        hist.SetDirectory(0)
        histograms[name] = hist

    return histograms


def create_pdf_histograms(count_histograms: dict[str, TH1F], category: str) -> dict[str, TH1F]:
    """Create PDFs from raw counts. Never use these as the resumable state."""
    pdf_histograms: dict[str, TH1F] = {}

    for name, count_hist in count_histograms.items():
        pdf_hist = count_hist.Clone(f"{category.lower()}_{name}_pdf")
        pdf_hist.SetDirectory(0)
        pdf_hist.SetTitle(count_hist.GetTitle().replace("Counts", "Probability Density"))

        integral = pdf_hist.Integral()
        if integral > 0:
            pdf_hist.Scale(1.0 / integral, "width")

        pdf_histograms[name] = pdf_hist

    return pdf_histograms


def check_normalization(signal_pdfs: dict[str, TH1F], background_pdfs: dict[str, TH1F]) -> None:
    print("\nChecking histogram normalization:")

    for category, histograms in [("Signal", signal_pdfs),("Background", background_pdfs)]:
        for name, hist in histograms.items():
            integral = hist.Integral("width")
            print(f"{category:10s} | {name:40s} | Integral = {integral:.6f}")



def _write_count_histograms(outfile, signal_histograms, background_histograms) -> None:
    signal_dir = outfile.mkdir("Signal")
    background_dir = outfile.mkdir("Background")

    signal_counts_dir = signal_dir.mkdir("Counts")
    background_counts_dir = background_dir.mkdir("Counts")

    signal_counts_dir.cd()
    for name, hist in signal_histograms.items():
        hist.Write(name)

    background_counts_dir.cd()
    for name, hist in background_histograms.items():
        hist.Write(name)


def write_final_histograms(signal_histograms: dict[str, TH1F], background_histograms: dict[str, TH1F]) -> None:
    """Write counts and derived PDFs to the final output, using atomic replacement."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    signal_pdfs = create_pdf_histograms(signal_histograms, "signal")
    background_pdfs = create_pdf_histograms(background_histograms, "background")
    check_normalization(signal_pdfs, background_pdfs)

    outfile = ROOT.TFile.Open(str(FINAL_OUTPUT_FILE), "RECREATE")

    try:
        _write_count_histograms(outfile, signal_histograms, background_histograms)

        signal_pdf_dir = outfile.GetDirectory("Signal").mkdir("PDFs")
        background_pdf_dir = outfile.GetDirectory("Background").mkdir("PDFs")

        signal_pdf_dir.cd()
        for name, hist in signal_pdfs.items():
            hist.Write(name)

        background_pdf_dir.cd()
        for name, hist in background_pdfs.items():
            hist.Write(name)

        outfile.cd()
        outfile.Write()
    finally:
        outfile.Close()


def main():
    # read the checkpoint file
    checkpoint_file = ROOT.TFile.Open(str(OUTPUT_DIR / "eID_MC_checkpoint.root"), "READ")
    hist_names = ["pT_hist", "EoverP_hist", "normalized_electron_likelihood_hist", "isolation_energy_ratio_hist"]
    signal_histograms = read_histogram_set(checkpoint_file, "Signal", hist_names)
    background_histograms = read_histogram_set(checkpoint_file, "Background", hist_names)

    write_final_histograms(signal_histograms, background_histograms)

if __name__ == "__main__":
    main()

