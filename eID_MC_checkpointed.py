import json
import os
from pathlib import Path

import ROOT  # pyright: ignore[reportMissingImports]
from ROOT import TH1F  # pyright: ignore[reportMissingImports]
from ROOT.Math import PxPyPzEVector  # pyright: ignore[reportMissingImports]

ROOT.TH1.AddDirectory(False)

from Libs.particle_calculations import calculate_elctron_likelihood, calculate_isolation_energy
from Libs.file_reader_checkpointed import get_file_paths, read_single_file  

    
# ---------------------------------------------------------------------------
# Run configuration
# ---------------------------------------------------------------------------
NUMBER_OF_FILES = 100  # In streaming mode this means 100 run numbers per Q2 bin.
STREAM_INPUT = True

# Save a checkpoint after this many events within each input file.
CHECKPOINT_EVERY_EVENTS = 1_000

OUTPUT_DIR = Path("Output_Root")
CHECKPOINT_FILE = OUTPUT_DIR / "eID_MC_checkpoint.root"
FINAL_OUTPUT_FILE = OUTPUT_DIR / "eID_MC_hists.root"


# ---------------------------------------------------------------------------
# Histogram creation and ROOT I/O
# ---------------------------------------------------------------------------

def _make_histogram_set(prefix: str) -> dict[str, TH1F]:
    nbins = 500
    histograms = {
        "pT_hist": TH1F(
            f"{prefix}_pT_hist",
            "Transverse Momentum Distribution; pT (GeV/c); Counts",
            nbins,
            0,
            20,
        ),
        "EoverP_hist": TH1F(
            f"{prefix}_EoverP_hist",
            "E/p Distribution; E/p; Counts",
            nbins,
            0,
            3,
        ),
        "normalized_electron_likelihood_hist": TH1F(
            f"{prefix}_normalized_electron_likelihood_hist",
            "Normalized Electron Likelihood Distribution; Normalized Electron Likelihood; Counts",
            nbins,
            0,
            1.1,
        ),
        "isolation_energy_ratio_hist": TH1F(
            f"{prefix}_isolation_energy_ratio_hist",
            "Isolation Energy Ratio Distribution; Isolation Energy Ratio; Counts",
            nbins,
            0,
            1.1,
        ),
    }

    # Explicit sum-of-squares storage is useful if event weights are introduced later.
    for hist in histograms.values():
        hist.Sumw2()

    return histograms


def create_count_histograms() -> tuple[dict[str, TH1F], dict[str, TH1F]]:
    signal = _make_histogram_set("signal")
    background = _make_histogram_set("background")
    return signal, background


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


def _read_histogram_set(root_file, directory: str, names: list[str]) -> dict[str, TH1F]:
    histograms: dict[str, TH1F] = {}

    for name in names:
        stored_hist = root_file.Get(f"{directory}/Counts/{name}")
        if not stored_hist:
            raise RuntimeError(
                f"Checkpoint is missing histogram {directory}/Counts/{name}"
            )

        # Clone so the histograms remain valid after closing the ROOT input file.
        hist = stored_hist.Clone(f"{directory.lower()}_{name}")
        hist.SetDirectory(0)
        histograms[name] = hist

    return histograms


def save_checkpoint(
    signal_histograms: dict[str, TH1F],
    background_histograms: dict[str, TH1F],
    metadata: dict,
) -> None:
    """Atomically replace the checkpoint with a complete cumulative snapshot."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    temporary_file = CHECKPOINT_FILE.with_name(CHECKPOINT_FILE.stem + ".tmp.root")

    outfile = ROOT.TFile.Open(str(temporary_file), "RECREATE")
    if not outfile or outfile.IsZombie():
        raise OSError(f"Could not create checkpoint file: {temporary_file}")

    try:
        _write_count_histograms(outfile, signal_histograms, background_histograms)
        outfile.cd()
        ROOT.TObjString(json.dumps(metadata)).Write("CheckpointMetadata")
        outfile.Write()
    finally:
        outfile.Close()

    # Both paths are in the same directory/filesystem. The previous valid checkpoint
    # remains in place until the new ROOT file has been completely written and closed.
    os.replace(temporary_file, CHECKPOINT_FILE)


def load_checkpoint(expected_input_files: list[str]):
    infile = ROOT.TFile.Open(str(CHECKPOINT_FILE), "READ")
    if not infile or infile.IsZombie():
        raise OSError(f"Could not open checkpoint file: {CHECKPOINT_FILE}")

    try:
        metadata_object = infile.Get("CheckpointMetadata")
        if not metadata_object:
            raise RuntimeError("Checkpoint does not contain CheckpointMetadata")

        metadata = json.loads(str(metadata_object.GetString().Data()))

        # Never resume against a different file list/order: doing so could silently
        # skip or double-count input data.
        if metadata.get("input_files") != expected_input_files:
            raise RuntimeError(
                "The input file list differs from the one stored in the checkpoint. "
                "Use the original configuration to resume, or remove the checkpoint "
                "to intentionally start a new run."
            )

        histogram_names = [
            "pT_hist",
            "EoverP_hist",
            "normalized_electron_likelihood_hist",
            "isolation_energy_ratio_hist",
        ]
        signal_histograms = _read_histogram_set(infile, "Signal", histogram_names)
        background_histograms = _read_histogram_set(infile, "Background", histogram_names)

    finally:
        infile.Close()

    return signal_histograms, background_histograms, metadata


def create_pdf_histograms(count_histograms: dict[str, TH1F], category: str) -> dict[str, TH1F]:
    """Create PDFs from raw counts. Never use these as the resumable state."""
    pdf_histograms: dict[str, TH1F] = {}

    for name, count_hist in count_histograms.items():
        pdf_hist = count_hist.Clone(f"{category.lower()}_{name}_pdf")
        pdf_hist.SetDirectory(0)
        pdf_hist.SetTitle(
            count_hist.GetTitle().replace("; Counts", "; Probability Density")
        )

        integral = pdf_hist.Integral()
        if integral > 0:
            pdf_hist.Scale(1.0 / integral, "width")

        pdf_histograms[name] = pdf_hist

    return pdf_histograms


def check_normalization(signal_pdfs: dict[str, TH1F], background_pdfs: dict[str, TH1F]) -> None:
    print("\nChecking histogram normalization:")

    for category, histograms in [
        ("Signal", signal_pdfs),
        ("Background", background_pdfs),
    ]:
        for name, hist in histograms.items():
            integral = hist.Integral("width")
            print(f"{category:10s} | {name:40s} | Integral = {integral:.6f}")


def write_final_histograms(
    signal_histograms: dict[str, TH1F],
    background_histograms: dict[str, TH1F],
    metadata: dict,
) -> None:
    """Write counts and derived PDFs to the final output, using atomic replacement."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    temporary_file = FINAL_OUTPUT_FILE.with_name(FINAL_OUTPUT_FILE.stem + ".tmp.root")

    signal_pdfs = create_pdf_histograms(signal_histograms, "signal")
    background_pdfs = create_pdf_histograms(background_histograms, "background")
    check_normalization(signal_pdfs, background_pdfs)

    outfile = ROOT.TFile.Open(str(temporary_file), "RECREATE")
    if not outfile or outfile.IsZombie():
        raise OSError(f"Could not create final ROOT file: {temporary_file}")

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
        ROOT.TObjString(json.dumps(metadata)).Write("RunMetadata")
        outfile.Write()
    finally:
        outfile.Close()

    os.replace(temporary_file, FINAL_OUTPUT_FILE)


# ---------------------------------------------------------------------------
# Event processing -- physics selections/calculations retained from your code
# ---------------------------------------------------------------------------

def process_event(
    frame,
    signal_histograms: dict[str, TH1F],
    background_histograms: dict[str, TH1F],
    counters: dict[str, int],
) -> None:
    reco_particles = frame.get("ReconstructedParticles")
    scattered_electrons_truth = frame.get("ScatteredElectronsTruth")

    counters["num_reco_particles"] += len(reco_particles)
    num_scattered_electrons = len(scattered_electrons_truth)

    if num_scattered_electrons == 0:
        return

    # Your original code uses the first truth scattered electron for matching.
    scattered_electron_id = scattered_electrons_truth[0].getObjectID().index

    for reco_particle in reco_particles:
        charge = reco_particle.getCharge()
        number_of_clusters = reco_particle.getClusters().size()
        number_of_tracks = reco_particle.getTracks().size()

        if number_of_clusters == 0 or number_of_tracks == 0:
            continue

        if charge != -1:
            continue

        momentum = reco_particle.getMomentum()
        energy = reco_particle.getEnergy()
        mom4 = PxPyPzEVector(momentum.x, momentum.y, momentum.z, energy)

        cluster = reco_particle.getClusters()[0]
        cluster_energy = cluster.getEnergy()
        cluster_position = cluster.getPosition()

        EoverP = cluster_energy / mom4.P()
        isolation_energy = calculate_isolation_energy(cluster_position, reco_particles)
        isolation_energy_ratio = cluster_energy / isolation_energy

        pids = reco_particle.getParticleIDs()
        number_of_pid = pids.size()

        pid_flag = True
        if number_of_pid > 0:
            normalized_electron_likelihood = calculate_elctron_likelihood(reco_particle)
        else:
            pid_flag = False
            normalized_electron_likelihood = None

        reco_id = reco_particle.getObjectID().index

        if scattered_electron_id == reco_id:
            counters["num_scattered_electrons_total"] += 1
            signal_histograms["pT_hist"].Fill(mom4.Pt())
            signal_histograms["EoverP_hist"].Fill(EoverP)
            signal_histograms["isolation_energy_ratio_hist"].Fill(isolation_energy_ratio)

            if pid_flag:
                signal_histograms["normalized_electron_likelihood_hist"].Fill(
                    normalized_electron_likelihood
                )
        else:
            counters["num_background_particles_total"] += 1
            background_histograms["pT_hist"].Fill(mom4.Pt())
            background_histograms["EoverP_hist"].Fill(EoverP)
            background_histograms["isolation_energy_ratio_hist"].Fill(isolation_energy_ratio)

            if pid_flag:
                background_histograms["normalized_electron_likelihood_hist"].Fill(
                    normalized_electron_likelihood
                )


# ---------------------------------------------------------------------------
# Main processing and restart logic
# ---------------------------------------------------------------------------

def main() -> None:
    input_files = get_file_paths(NUMBER_OF_FILES, stream=STREAM_INPUT, debug=False)
    print(f"Input files: {len(input_files)}")

    counters = {
        "num_reco_particles": 0,
        "num_scattered_electrons_total": 0,
        "num_background_particles_total": 0,
        "events_processed": 0,
    }

    if CHECKPOINT_FILE.exists():
        print(f"Loading checkpoint: {CHECKPOINT_FILE}")
        signal_histograms, background_histograms, metadata = load_checkpoint(input_files)
        counters.update(metadata["counters"])
        next_file_index = int(metadata["next_file_index"])
        next_event_index = int(metadata["next_event_index"])

        completed_files = metadata.get("completed_files", [])
        print(f"Resuming: {len(completed_files)}/{len(input_files)} files completed")
        if next_file_index < len(input_files):
            print(
                f"Current file: {input_files[next_file_index]}\n"
                f"Next event index in this file: {next_event_index}"
            )
    else:
        print("No checkpoint found; starting a new run.")
        signal_histograms, background_histograms = create_count_histograms()
        next_file_index = 0
        next_event_index = 0
        completed_files = []
        metadata = {}

    for file_index in range(next_file_index, len(input_files)):
        file_path = input_files[file_index]
        resume_event_index = next_event_index if file_index == next_file_index else 0

        print(f"\n[{file_index + 1}/{len(input_files)}] Opening input file:")
        print(file_path)

        # Keep reader referenced until all frames from this file have been processed.
        reader, frames = read_single_file(file_path)

        for event_index, frame in enumerate(frames):
            if event_index < resume_event_index:
                continue

            process_event(frame, signal_histograms, background_histograms, counters)
            counters["events_processed"] += 1

            if (event_index + 1) % 100 == 0:
                print(
                    f"File event {event_index + 1}; "
                    f"total processed events: {counters['events_processed']}"
                )

            if (event_index + 1) % CHECKPOINT_EVERY_EVENTS == 0:
                metadata = {
                    "input_files": input_files,
                    "next_file_index": file_index,
                    "next_event_index": event_index + 1,
                    "completed_files": completed_files,
                    "current_file": file_path,
                    "counters": counters.copy(),
                    "checkpoint_every_events": CHECKPOINT_EVERY_EVENTS,
                    "completed": False,
                }
                save_checkpoint(signal_histograms, background_histograms, metadata)
                print(
                    f"Checkpoint saved after event {event_index + 1} of this file "
                    f"({counters['events_processed']} events total)."
                )

        # This file has been fully read. Mark it complete in the same checkpoint
        # that stores the corresponding cumulative histogram contents.
        completed_files = input_files[: file_index + 1]
        metadata = {
            "input_files": input_files,
            "next_file_index": file_index + 1,
            "next_event_index": 0,
            "completed_files": completed_files,
            "current_file": None,
            "counters": counters.copy(),
            "checkpoint_every_events": CHECKPOINT_EVERY_EVENTS,
            "completed": False,
        }
        save_checkpoint(signal_histograms, background_histograms, metadata)
        print(f"Finished input file {file_index + 1}/{len(input_files)}.")

        # Avoid retaining the previous file's event collection longer than needed.
        del frames
        del reader

        next_event_index = 0

    metadata = {
        "input_files": input_files,
        "next_file_index": len(input_files),
        "next_event_index": 0,
        "completed_files": input_files,
        "current_file": None,
        "counters": counters.copy(),
        "checkpoint_every_events": CHECKPOINT_EVERY_EVENTS,
        "completed": True,
    }
    save_checkpoint(signal_histograms, background_histograms, metadata)

    print("\nAnalysis summary:")
    print("Total number of events processed: ", counters["events_processed"])
    print("Total number of reconstructed particles: ", counters["num_reco_particles"])
    print(
        "Total number of scattered electrons matched to reconstructed particles: ",
        counters["num_scattered_electrons_total"],
    )
    print(
        "Total number of background particles in events with scattered electron: ",
        counters["num_background_particles_total"],
    )

    write_final_histograms(signal_histograms, background_histograms, metadata)
    print(f"\nFinal histograms written to: {FINAL_OUTPUT_FILE}")
    print(f"Restart checkpoint retained at: {CHECKPOINT_FILE}")


if __name__ == "__main__":
    main()
