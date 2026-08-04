Code for My Marine Ecosystem Managment MSc Thesis - Illuminating Reef architecture

master_script.R
main pipeline script - prepares DEM, computes direct and ambient shading, pulls Kd sources, and water depth, and for each logger and hour predicts PAR - and creates validation sheets against the measured values

big_shadow_maps.R
For each hourly timestep of chosen sampling event, raytraces the shadow over the reef DEM to output PNGs per timestep showing the shadow and light regimes

shadow_diagnostic.R
picks one logger and traces every component of pred par across the day (needs to have full pipeline run)


kdpar_full.R
pulls GBR spectral Kd data and collapses it to Kd lookup tables

