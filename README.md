Scripts (in progress) to perform sensor and source-space decoding analyses on Fieldtrip-processed MEG data (using a Support Vector Machine for two-class problems), using different spatiotemporal feature selection approaches and cross-validation schemes. 
Work in progress: scripts for performing Representational Similarity Analysis of spatiotemporally resolved MEG data. Reliant on [Fieldtrip](http://www.fieldtriptoolbox.org/) for preprocessing, plotting & templates and [LibLinear](https://www.csie.ntu.edu.tw/~cjlin/liblinear/) for SVM decoding.

# Structure

## [Prepare](https://github.com/dianadima/mvpa-for-meg/wiki/Prepare)
Scripts for preprocessing sensor-level MEG data, source-reconstructing MEG data with a beamformer, and preparing data for decoding (creating pseudo-trials, whitening, getting spatial clustering information for spatially-resolved MVPA).

## [Decode](https://github.com/dianadima/mvpa-for-meg/wiki/Decode)
Scripts for decoding MEG data across time (time-resolved and temporal generalization) and space (searchlight, region-of-interest, source/channel selection, whole-brain), with kfold or hold-out cross-validation.
The main scripts use the LibLinear SVM. Scripts based on LibSVM and Matlab SVM are included in subfolders but have not been updated/tested...

## Plot
Scripts for plotting decoding results over time and in space (sensor/source) and as movies (over time/rotating brain).

## Stats
Non-parametric significance testing with different methods of correcting for multiple comparisons.

## RSA
Representational Similarity Analysis scripts (work in progress): create dissimilarity matrices, correlate them with models, calculate a noise ceiling and randomize the correlations for significance testing. 
