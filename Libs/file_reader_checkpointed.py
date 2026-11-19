from podio import Frame  # pyright: ignore[reportMissingImports]
from podio.reading import get_reader  # pyright: ignore[reportMissingImports]


def get_file_paths(n_files: int, stream: bool = False, debug: bool = False) -> list[str]:
    """Return the input file paths in the same order as the original reader."""
    if n_files < 1:
        raise ValueError("Number of files must be at least 1.")

    paths: list[str] = []

    if stream:
        for i in range(n_files):
            paths.extend(
                [
                    f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_1to10/pythia8.316-1.0_NC_noRad_ep_10x100_q2_1to10_run000.{i:04d}.eicrecon.edm4eic.root",
                    f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_10to100/pythia8.316-1.0_NC_noRad_ep_10x100_q2_10to100_run000.{i:04d}.eicrecon.edm4eic.root",
                    f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_100to1000/pythia8.316-1.0_NC_noRad_ep_10x100_q2_100to1000_run000.{i:04d}.eicrecon.edm4eic.root",
                    f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_1000toINF/pythia8.316-1.0_NC_noRad_ep_10x100_q2_1000toINF_run000.{i:04d}.eicrecon.edm4eic.root",
                ]
            )
    else:
        if n_files > 3:
            raise ValueError(
                "Number of local files must not exceed 3. "
                "I have not downloaded more than 3 files yet."
            )

        for i in range(n_files):
            paths.append(
                f"Data/pythia8.316-1.0_NC_noRad_ep_10x100_q2_10to100_run001.{i:04d}.eicrecon.edm4eic.root"
            )

    if debug:
        print(f"Found {len(paths)} input files")

    return paths


def read_single_file(file_path: str):
    """Open one EDM4eic file and return both reader and its event collection.

    Keep the returned reader alive while iterating over the returned events.
    """
    reader = get_reader([file_path])
    frames = reader.get("events")
    return reader, frames


def file_reader(n_files: int, stream: bool = False, debug: bool = False) -> list[Frame]:
    """Compatibility wrapper matching the original all-files-at-once interface."""
    paths = get_file_paths(n_files=n_files, stream=stream, debug=debug)
    reader = get_reader(paths)
    frames = reader.get("events")
    return frames
