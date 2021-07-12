This repository contains code and data for the paper _A Random-Walk-Based Epidemiological Model_ by Andrew Chu, Greg Huber, Aaron McGeever, Boris Veytsman and David Yllanes

## Directories and files ##

* [Manifest.toml](Manifest) - used to set up Julia environment
* [Project.toml](Project) - used to set up Julia environment
* [data](data) - data files. ".jld" files are compatible with Julia. Additional data is provided in this [Google Drive folder](https://drive.google.com/drive/folders/1FtTw2a9yRfbeCAfFGvMAgRWLOo2wszKb?usp=sharing) due to Github file size limits.
* [R0-prediction](R0-prediction) - code for generation of R0 prediction figure
* [Figure 1c](Figure_1c.ipynb) - code for generating figure 1c
* [Figure 3](Figure_3.ipynb) - code for generating figure 3
* [Figure 4](Figure_4.ipynb) - code for generating figure 4
* [Figure 5](Figure_5.ipynb) - code for generating figure 5. Uses data from "Figure_5_data.jld" in data folder
* Figure 6 - Uses data located in the [Google Drive folder](https://drive.google.com/drive/folders/1FtTw2a9yRfbeCAfFGvMAgRWLOo2wszKb?usp=sharing).
* Figure 7 - Uses data located in the [Google Drive folder](https://drive.google.com/drive/folders/1FtTw2a9yRfbeCAfFGvMAgRWLOo2wszKb?usp=sharing).
* [Table 1](Table_1.ipynb) - code for generating table 1
* [Figure S2](Figure_S2.ipynb) - code for generating figure S2. Uses data from files beginning with "rg_ar_" in data folder
* [Figures S3 - S5](Figures_S3-5.ipynb) - code for generating figures S3 - S5. Uses data located in the [Google Drive folder](https://drive.google.com/drive/folders/1FtTw2a9yRfbeCAfFGvMAgRWLOo2wszKb?usp=sharing).
* [Iso-R0 Estimation](isoR0Estimation.ipynb) - code for generating data for iso-R0 lines in figure 4
* [Phase Diagram](phasediagram.ipynb) - code for generating data for phase diagram in figure 4
* [Radius of Gyration and Attack Rate Calculation](rofgyr_attrate_calculation.ipynb) - code for calculating radius of gyration and attack rate data in figure 6.
* [Cluster Simulation](clusterSimulation.jl) - code for simulating outbreaks and calculating the radius of gyration, cluster mass, and cluster surface points in figures 6 and 7.
* [utils.jl](utils.jl) - helper functions for [Cluster Simulation](clusterSimulation.jl).
