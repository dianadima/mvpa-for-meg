function [ rand_accuracy ] = rand_searchlight_kfold( data, labels, cluster_idx, num_iterations, varargin )
% Inputs: data, labels, cluster_idx (neighbourhood structure or source indices obtained using get_sensor_info or get_source_info).
% Optional: 
%          'channels', channel set set (string or cell array of strings; default: 'MEG'). 
%          'decoding_window' (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%           other possible name-value pairs: SVM settings, svm evaluation metrics (see Documentation).
% Output: time and space-resolved accuracy and F1-score (channels x time).
% Performs time and space-resolved SVM decoding of MEG data, using stratified k-fold cross-validation and LibLinear for each time window and sensor cluster.
% Uses sensor neighbours defined using FT template configuration (or custom structure containing source indices).
%
% DC Dima 2018 (diana.c.dima@gmail.com)

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

parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','decoding_window', 'time'}); %converted struct will be fed into decoding function
clear p;

%channel indices for each searchlight iteration
if isstruct(cluster_idx) %the sensor-space case
    chan_idx = arrayfun(@(i) find(ismember({cluster_idx.label},[cluster_idx(i).label; cluster_idx(i).neighblabel])), 1:length(cluster_idx), 'UniformOutput', false); %store all searchlight idx in a cell array
else
    chan_idx = cluster_idx; %the source-space case
end;

if isempty(dec_args.cv_indices)
    cv = cvpartition(labels, 'kfold', 5);
    dec_args.cv_indices = cv;
end;

fprintf('\nRunning %d searchlights...\n', length(chan_idx)); 

data_svm = arrayfun(@(i) data(chan_idx{i}, :,:), 1:length(chan_idx), 'UniformOutput', false); %channel and time selection
accuracy = arrayfun(@(i) time_resolved_kfold_rand(data_svm{i}, labels, num_iterations, dec_args, svm_par), 1:length(data_svm), 'UniformOutput', false);
rand_accuracy = cell2mat([accuracy]);

end

