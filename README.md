## Bayesian method for DIS Electron ID for ePIC experiment

Repo used as development codebase for project on scattered electron ID using Bayesian PID methods.
Assumes running in current EIC shell environemnt and uses edm4eic data_model.

### Analysis overview

- **MC loop (`eID_MC.py`):** Matches reconstructed particles to the truth scattered electron and builds signal and background distributions for transverse momentum, E_cal-to-P_track ratio, isolation energy threshold for Calorimeter cluster, and PID likelihood for PID detectors.
- **Reconstructed loop (`eID_Rec.py`):** Uses those MC distributions in a Bayesian calculation to assign each selected reconstructed particle a scattered-electron probability.
- **Plotter (`plot_podio.py`):** Saves a lot of histrograms about different podio collections used. Usually for my initial developement and just messing around 
