# The DESI Bright Galaxy Survey Early Data Analysis

`getting_started.ipynb`: getting started on NERSC

`desi_plots.ipynb`: some basic DESI plots and hopefully some helpful info

In each folder of SV3, DA0.2:
  - `bright_cat.ipynb`: collection of plots looking into the properties of each sample
  - `k+e_vol_lim.ipynb`: using k+e corrections to define volume-limited samples used for clustering.

`clustering_plots.ipynb`: collection of clustering plots for both SV3 and DA0.2

`plotting.py`: makes histograms from stepwise into nice curves
`vol_lim_samples.py`: functions for defining volume-limited samples 

copy of my version of `xirunpc.py` - original found at https://github.com/desihub/LSS/blob/main/scripts/xirunpc.py - I have made additions to allow cutting to volume-limited samples and splitting into colour (needs to be placed in the appropriate location in LSS directory to run)

copy of my version of `mkCat_main.py` - original found at https://github.com/desihub/LSS/blob/main/scripts/main/mkCat_main.py - edited by Sam and myself to work with the version of k-correction code at the time (propbably now out of date?)