function [ accuracy, Fscore] = searchlight_kfold_libsvm( data, labels, cluster_idx, varargin )
% Performs time-resolved SVM decoding of MEG data, using stratified k-fold cross-validation on each time window, and LibSVM SVM implementation.
% Inputs: data, labels. Data format: channels/sources x time x trials
%         cluster_idx, structure obtained using get_sensor_info or get_source_info - for channel/source selection. You can also just provide numerical indices, in which case you don't need the structure.
% Optional: 
%          'channels', channel set (string or cell array of strings; default: 'MEG').
%          'decoding_window' (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%
%           Below are other name-value pairs that control the SVM settings, with defaults:
%           
%          boxconstraint = 1; --> C-parameter: Note, we don't have any options for optimizing this, need to write it separately if needed
%          kfold = 5; --> number of folds for k-fold-cross-validation
%          cv_indices = []; --> supply cross-validation indices (e.g. in cvpartition format)
%          standardize = true; --> standardize features using mean and SD of training set (recommended)
%          weights = false; --> calculate weights (by retraining model on whole dataset)
%          kernel = 0; --> only applies to libSVM: 	0: linear; 1: polynomial; 2:radial basis function; 3:sigmoid
%
% Outputs: structure containing classification performance metrics (for each timepoint and cross-validation round).
%
% DC Dima 2018 (diana.c.dima@gmail.com)

addpath(strrep(mfilename('fullpath'),fullfile('libsvm','searchlight_kfold_libsvm'),''));

%parse inputs
dec_args = decoding_args;
svm_par = svm_args;
list = [fieldnames(dec_args); fieldnames(svm_par)];
p = inputParser;

for i = 1:length(properties(decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
for ii = i+1:length(properties(dec_args))+length(properties(svm_args))
    addParameter(p, list{ii}, svm_par.(list{ii}));
end;

addParameter(p,'time',[]);


parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','decoding_window', 'time','mnn','pseudo'}); %converted struct will be fed into decoding function
clear p;

%channel indices for each searchlight iteration
if isstruct(cluster_idx) %the sensor-space case
    chan_idx = arrayfun(@(i) find(ismember({cluster_idx.label},[cluster_idx(i).label; cluster_idx(i).neighblabel])), 1:length(cluster_idx), 'UniformOutput', false); %store all searchlight idx in a cell array
else
    chan_idx = cluster_idx; %the source-space case
end;

%create time axis
if ~isempty(dec_args.time)
    time = dec_args.time;
elseif isfield(cluster_idx, 'time')
    time = cluster_idx.time;
else
    time = 1:size(data,2);
end;

%time limits for decoding window
if ~isempty(dec_args.decoding_window)
    if ~isempty(find(round(time,3)==dec_args.decoding_window(1),1))
        lims(1) = find(round(time,3)==dec_args.decoding_window(1));
    else
        fprintf('\nWarning: starting timepoint not found, starting from beginning of data...\n');
        lims(1) = 1;
    end;
    if ~isempty(find(round(time,3)==dec_args.decoding_window(2),1))
        lims(2) = find(round(time,3)==dec_args.decoding_window(2));
    else
        fprintf('\nWarning: end timepoint not found, decoding until end of data...\n');
        lims(2) = size(data,2);
    end;
        
else
    lims = [1 size(data,2)];
end;

if size(data,2)<lims(2)
    lims(2) = size(data,2);
end;

accuracy = zeros(length(chan_idx), length(lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1));
Fscore = zeros(length(chan_idx), length(lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1));
fprintf('\nRunning searchlight '); 

for c = 1:length(chan_idx)
    
    fprintf('%d out of %d', c, length(chan_idx));
    data_svm = arrayfun(@(i) reshape(data(chan_idx{c}, i:i+dec_args.window_length-1,:), length(chan_idx{c})*dec_args.window_length, size(data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %channel and time selection
    results = arrayfun(@(i) svm_decode_kfold_libsvm(data_svm{i},labels,svm_par,'AUC', false), 1:length(data_svm));
    accuracy(c,:) = cell2mat({results.Accuracy});
    Fscore(c,:) = cell2mat({results.WeightedFscore});
    fprintf((repmat('\b',1,numel([num2str(c) num2str(length(chan_idx))])+8)));
    clear data_svm;
    
end;

end

