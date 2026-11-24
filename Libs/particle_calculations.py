#import edm4hep                                  # pyright: ignore[reportMissingImports]
from edm4hep import ReconstructedParticle       # pyright: ignore[reportMissingImports]
import math               
from ROOT.Math import PxPyPzEVector, XYZVector             # pyright: ignore[reportMissingImports]

PDG_CODE_ELECTRON = 11
PDG_CODE_PION = 211
PDG_CODE_KAON = 321
PDG_CODE_PROTON = 2212


def calculate_elctron_likelihood(reco_particle: ReconstructedParticle) -> list[float]:
    pids = reco_particle.getParticleIDs()
    max_likelihoods = {"electron": 0.0, "pion": 0.0, "kaon": 0.0, "proton": 0.0}
    normalized_electron_likelihoods = []

    for i, pid in enumerate(pids):
        assigned_pdg = abs(pid.getPDG())
        likelihood = pid.getLikelihood()

        if i%4 == 0 and i > 0:  # Clear the max_likelihoods dictionary every 4 iterations to avoid accumulating likelihoods from different PID objects. This is based on the assumption that each ReconstructedParticle has 4 PID objects corresponding to electron, pion, kaon, and proton.
            print(f"Likelihoods so far: {max_likelihoods} for pid index {i}, compute electron likelihood and then clear for next PID object.")

            total_likelihood = sum(max_likelihoods.values())
            normalized_electron_likelihood = max_likelihoods["electron"] / total_likelihood if total_likelihood > 0 else -1.0
            print(f"Normalized electron likelihood: {normalized_electron_likelihood}, total likelihood: {total_likelihood}")
            normalized_electron_likelihoods.append(normalized_electron_likelihood)

            max_likelihoods = {"electron": 0.0, "pion": 0.0, "kaon": 0.0, "proton": 0.0}


        if assigned_pdg == PDG_CODE_ELECTRON:  # Electron
            max_likelihoods["electron"] = max(max_likelihoods["electron"], likelihood)
        elif assigned_pdg == PDG_CODE_PION:  # Pion
            max_likelihoods["pion"] = max(max_likelihoods["pion"], likelihood)
        elif assigned_pdg == PDG_CODE_KAON:  # Kaon
            max_likelihoods["kaon"] = max(max_likelihoods["kaon"], likelihood)
        elif assigned_pdg == PDG_CODE_PROTON:  # Proton
            max_likelihoods["proton"] = max(max_likelihoods["proton"], likelihood)
        else:
            print(f"Warning: Unrecognized PDG code {assigned_pdg} encountered. Ignoring this PID.")
        
    
    return normalized_electron_likelihoods
    


def calculate_isolation_energy(cluster_position, reconstructed_particles, isolation_radius=0.8) -> float:
    isolation_energy = 0.0

    lead_pos = XYZVector(cluster_position.x, cluster_position.y, cluster_position.z)
    lead_eta = lead_pos.Eta()
    lead_phi = lead_pos.Phi()

    for reco_particle in reconstructed_particles:
        for cluster in reco_particle.getClusters():
            cluster_pos = XYZVector(cluster.getPosition().x, cluster.getPosition().y, cluster.getPosition().z)
            cluster_eta = cluster_pos.Eta()
            cluster_phi = cluster_pos.Phi()

            cone_radius = math.sqrt((lead_eta - cluster_eta) ** 2 + (lead_phi - cluster_phi) ** 2)

            if cone_radius < isolation_radius:
                isolation_energy += cluster.getEnergy()

    return isolation_energy