### Overview
This repository contains code and data to reproduce the results in Lesk & Mankin, 2026. 

The main aspects of the analysis are:
1) processing of raw data products (obtainable from URLs in Table S1) into derived data products, such as water-year daily precipitation Gini index
2) processing and assembling raw and derived data into the panels
3) running panel regressions (for Figures 1-3)
4) running the simple land-surface model (SEMB-H) (for Figure 4)
5) calculating projections (for Figure 5, ) 
6) plotting main and SI figures

Please find more detailed workflow overviews for each aspect below. From the data and code provided in this repository, you can reproduce the results from the raw data or from various intermediate steps. For example, you can simply reproduce the figures from the included regression and SEMB-H modelling result data, or rerun the regressions and SEMB-H models, or start from scratch and recompile the panels from the raw data.
### Workflow overview

Data processing, analysis, and plotting are coded up in Jupyter notebooks. The environment .yaml file used for this analysis is included in the repository to facilitate package compatibility. The regression analysis is in R and the required packages are enumerated at the top of the code.

The core intermediate data in the analysis is the compiled panels of water-year anomalies of GRACE TWS, daily precipitation Gini index, CPC temperature, GEWEX-SRB solar radiation. These are included in the /panels/ directory of this repository. The code to assemble these panels is in 'Lesk+Mankin_data-prep.ipynb'. To reproduce the daily precipitation Gini indices, use 'Lesk+Mankin_gini-data.ipynb'.

From this point, the analysis moves to R to run the panel models. To do so, use the code 'Lesk+Mankin_regressions.R'.




