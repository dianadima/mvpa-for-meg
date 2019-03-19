Functions for performing sensor and source-space decoding analyses on Fieldtrip-processed MEG data using a Support Vector Machine for two-class problems, implementing different spatiotemporal feature selection approaches and cross-validation schemes.  

Work in progress: functions for performing Representational Similarity Analysis of spatiotemporally resolved MEG data. 

Reliant on the [Fieldtrip](http://www.fieldtriptoolbox.org/) toolbox for preprocessing, source localization, plotting, and templates, and [LibLinear](https://www.csie.ntu.edu.tw/~cjlin/liblinear/) for SVM decoding.

# [Demo](https://github.com/dianadima/mvpa-for-meg/wiki/Demo)

An example dataset and script for running some decoding analyses are included in the [demo](https://github.com/dianadima/mvpa-for-meg/tree/master/demo) directory and described [here](https://github.com/dianadima/mvpa-for-meg/wiki/Demo).

# Structure

## [Prepare](https://github.com/dianadima/mvpa-for-meg/wiki/1.-Prepare)
Functions for preprocessing sensor-level MEG data, source-reconstructing MEG data with a beamformer, and preparing data for decoding (creating pseudo-trials, whitening, getting spatial clustering information for spatially-resolved MVPA).

## [Decode](https://github.com/dianadima/mvpa-for-meg/wiki/2.-Decode)
Functions for decoding MEG data across time (time-resolved and temporal generalization) and space (searchlight, region-of-interest, source/channel selection, whole-brain), with kfold or hold-out cross-validation.

The main functions use the LibLinear SVM. Functions based on LibSVM and Matlab SVM are included in subfolders.

## [Stats](https://github.com/dianadima/mvpa-for-meg/wiki/3.-Stats)
Non-parametric significance testing with different methods of correcting for multiple comparisons.

## [Plot](https://github.com/dianadima/mvpa-for-meg/wiki/4.-Plot)
Functions for plotting decoding results over time and in space (sensor/source space) and as movies (over time/on a rotating template brain).

## [RSA](https://github.com/dianadima/mvpa-for-meg/wiki/5.-RSA)
Functions for Representational Similarity Analysis (work in progress): create dissimilarity matrices, correlate them with models, calculate a noise ceiling and randomize the correlations for significance testing. 
