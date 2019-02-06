function [ results ] = time_resolved_kfold_libsvm( data, labels, varargin )
% Performs time-resolved SVM decoding of MEG data, using stratified k-fold cross-validation on each time window, and LibLinear SVM implementation.
% Inputs: data, labels
% Optional: 'sensor_idx', structure obtained using get_sensor_info - for channel selection ONLY. 
%                        You can also just provide numerical indices, in which case you don't need the structure.
%                        If you want to subselect features on source space data, you need to manually provide numerical indice (i.e. 'channels', [1:1000]).
%          'channels', channel set set (string or cell array of strings; default: 'MEG').
%          'decoding_window' (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%             
%           Below are other name-value pairs that control the SVM settings, with defaults:
%           
%          boxconstraint = 1; --> C-parameter: Note, we don't have any options for optimizing this, need to write it separately if needed
%          kfold = 5; --> number of folds for k-fold-cross-validation
%          cv_indices = []; --> supply cross-validation indices (e.g. in cvpartition format)
%          iterate_cv = 1; --> rounds of cross-validation to run
%          standardize = true; --> standardize features using mean and SD of training set (recommended)
%          weights = false; --> calculate weights (by retraining model on whole dataset)
%          kernel = 0; --> only applies to libSVM: 	0: linear; 1: polynomial; 2:radial basis function; 3:sigmoid
%          AUC = true (compute AUC: longer time)
%          plotROC = false (plot ROC curve).
%
% Outputs: structure containing classification performance metrics (for each timepoint and cross-validation round).
%
% DC Dima 2018 (diana.c.dima@gmail.com)

addpath(strrep(mfilename('fullpath'),fullfile('libsvm','time_resolved_kfold_libsvm'),'')); %add parent directory to path

%parse inputs
dec_args = args.decoding_args;
svm_par = args.svm_args;
list = [fieldnames(dec_args); fieldnames(svm_par)];
p = inputParser;

for i = 1:length(properties(args.decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
for ii = i+1:length(properties(args.dec_args))+length(properties(args.svm_args))
    addParameter(p, list{ii}, svm_par.(list{ii}));
end;
addParameter(p, 'AUC', true);
addParameter(p,'plotROC', false);
addParameter(p, 'sensor_idx', []);

parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','sensor_idx','decoding_window', 'time'}); %converted struct will be fed into decoding function
clear p;

%get channel indices and time axis. Numerical channel indices take priority
if ~iscell(dec_args.channels) && ~ischar(dec_args.channels)
    chan_idx = dec_args.channels;
end;
if ~isempty(dec_args.sensor_idx)
    if ~exist('chan_idx', 'var')
        chan_idx = 1:size(data,1); %initialize with entire array
        if ~strcmp (dec_args.channels, 'MEG') %if we need to subselect sensors
            chan = [];
            for i = 1:length(dec_args.channels)
                idx = cellfun('isempty',strfind({dec_args.sensor_idx.label},dec_args.channels{i}));
                chan = [chan chan_idx(~idx)]; %#ok<AGROW>
            end;
            chan_idx = chan;
        end;
    end;
else
    if ~exist('chan_idx', 'var')
        chan_idx = 1:size(data,1);
    end;
end;

%create time axis
if ~isempty(dec_args.time)
    time = dec_args.time;
elseif isfield(dec_args.sensor_idx, 'time')
    time = dec_args.sensor_idx.time;
else
    time = 1:size(data,2);
end;

if length(time)~=size(data,2)
    time = 1:size(data,2);
    fprintf('Warning: time axis does not match dataset size. Replacing with default time axis...');
end;


if ~isempty(dec_args.decoding_window)
    if ~isempty(find(round(time,3)==dec_args.decoding_window(1),1))
        lims(1) = find(round(time,3)==dec_args.decoding_window(1));
    else
        fprintf('Warning: starting timepoint not found, starting from beginning of data...\n');
        lims(1) = 1;
    end;
    if ~isempty(find(round(time,3)==dec_args.decoding_window(2),1))
        lims(2) = find(round(time,3)==dec_args.decoding_window(2));
    else
        fprintf('Warning: end timepoint not found, decoding until end of data...\n');
        lims(2) = size(data,2);
    end;
        
else
    lims = [1 size(data,2)];
end;

if size(data,2)<lims(2)
    lims(2) = size(data,2);
end;

%here we create a cell array containing classification data for each time window
data_svm = arrayfun(@(i) reshape(data(chan_idx,i:i+dec_args.window_length-1,:), length(chan_idx)*dec_args.window_length, size(data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %features are channelsxtime,observations are trials    
%and here we run the classifier
results = arrayfun(@(i) svm_decode_kfold_libsvm(data_svm{i},labels, svm_par, 'AUC', false), 1:length(data_svm));
  
end
