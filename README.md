# rsfc-fnirs
## These MATLAB codes reproduce the simulation results of the article "Correction of global physiology in resting-state functional near-infrared-spectroscopy" by Pradyumna Lanka, Heather Bortfeld & Theodore J. Huppert.

####These codes require the NIRS Brain AnalyzIR toolbox to run. The toolbox can be downloaded from  https://github.com/huppertt/nirs-toolbox. 
Then please add the toolbox to the MATLAB path

### Description of files

#### Dist_fig.m:
Generates figs 4,5,6 & 11 in the manuscript. This code compares several connectivity methods in simulated data with different statistical assumptions including, (i) Random data, (ii) Data with temporal autocorrelated noise, (iii) Data with both temporal autocorrelation and shared global signal, (iv) Data with temporal autocorrelation, shared global signal and motion artifacts.
Connectivity methods examined include robust and non-robust versions of (a) Pearson's correlation, (b) AR correlation (pre-whitened correlation), (c) Partial correlation (short sep channels), (d) Modified multivariate granger causality (MVGC) with the zeroth lag term.
The code generates the type-I error control curves,  (ii) receiver operating characteristic (ROC) curves as well as the expected vs actual null distribution plots.

#### Lag.m:
Generates figures 9 & 10 in the manuscript. This codes analyses the ability of MVGC, partial correlation and modified MVGC with zero lag to model lagged and zero lagged relationships between channels.

#### Parcorr.m:
Generates fig 7 in the manuscript. This code compares AR partial correlation methods with PCA filtering using short separation, long and both short & long channels. Both AR partial correlation with robust regression and without robust regression are also compared.

#### Grangers.m
Generates fig 8 in the manuscript. This code compares MVGC methods controlling for short separation, long and both short & long channels after PCA filtering. Both MVGC with robust regression and without robust regression are also compared.

#### parcorrfig.m
Generates fig 1 in the manuscript. This code visually compares Pearson's correlation, AR correlation & AR partial correlation with the ground truth in their effectiveness to correct for temporal autocorrelation and global signal-induced spatial covariance.

#### TSfig.m

Generates fig 2 in the manuscript. The code generates plots that show the raw signal intensity for the simulated data as well as the normalized covariance matrix for the channels. Rows show simulated data with just (A) temporal autocorrelation, (B) with both temporal autocorrelation and global systemic physiology, and (C) with temporal autocorrelation, global systemic physiology, and motion artifacts. It also generates the channel covariance matrix.


### Function folder:
grangercausality: A wrapper function that calls mymvgc
mymvgc: Function runs multivariate granger causality analysis. Depending on the parameters also runs modified MVGC.
partial_corr: Function runs partial correlation analysis for assessing zero-lag connectivity between channels
simData_connectivity_shortsep: Simulates Resting-state fNIRS data with a probe containing short separation channels. The simulated data mimics RS-fNIRS data and contains both global physiological noise and temporal autocorrelation with similar statistical properties as experimental RS-fNIRS data.
