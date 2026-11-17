from podio import Frame                         # pyright: ignore[reportMissingImports]
from podio.reading import get_reader            # pyright: ignore[reportMissingImports]

def file_reader(n_files: int, stream: bool = False, debug: bool = False) -> list[Frame]:

    string_list = []
    
    if stream:
        for i in range(n_files):
            file_path_1 = f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_1to10/pythia8.316-1.0_NC_noRad_ep_10x100_q2_1to10_run000.{i:04d}.eicrecon.edm4eic.root"
            file_path_2 = f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_10to100/pythia8.316-1.0_NC_noRad_ep_10x100_q2_10to100_run000.{i:04d}.eicrecon.edm4eic.root"
            file_path_3 = f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_100to1000/pythia8.316-1.0_NC_noRad_ep_10x100_q2_100to1000_run000.{i:04d}.eicrecon.edm4eic.root"
            file_path_4 = f"root://epicxrd1.sdcc.bnl.gov:1095//eic/EPIC//RECO/26.07.1/epic_craterlake/DIS/pythia8.316-1.0/NC/noRad/ep/10x100/q2_1000toINF/pythia8.316-1.0_NC_noRad_ep_10x100_q2_1000toINF_run000.{i:04d}.eicrecon.edm4eic.root"
            string_list.append(file_path_1)
            string_list.append(file_path_2)
            string_list.append(file_path_3)
            string_list.append(file_path_4)

    else:
        if n_files < 1:
            raise ValueError("Number of files must be at least 1.")
        elif n_files > 3:
            raise ValueError("Number of files must not exceed 3. I have not downloaded more than 3 files yet:)")

        for i in range(n_files):
            file_path = f"Data/pythia8.316-1.0_NC_noRad_ep_10x100_q2_10to100_run001.{i:04d}.eicrecon.edm4eic.root"
            string_list.append(file_path)

    reader = get_reader(string_list)
    frames = reader.get("events")

    return frames