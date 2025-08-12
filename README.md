Code to accompany manuscript "Sustained benefit of personalized closed-loop deep brain stimulation for major depressive disorder" by Kristin K. Sellers, Ankit Khambhati, Elissa J. Hamlat, Inhauck Choi, Daniela A. Astudillo Maya, Alexandra G. Tremblay-McGaw, Noah Stapper, Joncarmen Mergenthaler, Joline M. Fan, Edward F. Chang, and Andrew D. Krystal

# System Requirements:

- Matlab 2021b
- ClusterBorder function for Matlab from https://www.mathworks.com/matlabcentral/answers/769877-draw-lines-around-specific-regions-in-imagesc-plot
- Python 3.10.4
    - numpy==2.2.5
    - pandas==2.2.3
    - tqdm==4.67.1
    - matplotlib==3.10.1
    - seaborn==0.13.2

# Installation Guide

- Follow standard steps for Matlab installation -- this typically takes 1-2 hours
- Follow standard steps for Python intallation and packages listed above --  this typically takes 1-2 hours


# Instructions for Use

- KS_Presidio_SymptomAndDeviceMetricPlots.m does calculations and produces figures for manuscript Figure 1; Extended Figure 1; Extended Figure 2; and Extended Figure 3
    - Uses input data `PR01_ComprehensiveMADRS.mat`, `PR01_AmbulatoryRedcap.csv`, `PR01_Histogram_Hourly.csv`
    - Uses function `PresidioPatientData_PR01`

- KS_Presidio_RNS_Spectrograms_MADRS.m does calculations and produces figures for manuscript Figure 2B, C, D, E; Extended Figure 4; Extended Figure 5
    - Uses input data `PR01_ComprehensiveMADRS.mat`, `PR01_SelectMagnetData.mat`, and `PR01_LongitudinalSpectrograms.mat`
    - Uses function `PresidioPatientData_PR01`

- KS_Presidio_RNS_Biomarker.m does calculations and produces figures for manuscript Figure 2F; Extended Figure 4B
    - Uses input data `PR01_SelectData.mat`
    - Uses function `PresidioPatientData_PR01`

- AK_Presidio_MADRS_Decoder.ipynb is a Python-based Jupyter Notebook that does calculations and produces figures for manuscript Figure 2G.
    - Uses input data `PR01-LTFU-Neural_MADRS.csv`
