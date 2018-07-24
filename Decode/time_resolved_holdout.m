function [ results ] = time_resolved_holdout( train_data, train_labels, test_data, test_labels, varargin )
% Inputs: data, labels, cluster_idx (neighbourhood structure or source indices obtained using get_sensor_info or get_source_info).
% Data format: channels/sources x time x trials.
% Optional: 
%          'channels', channel set set (string or cell array of strings; default: 'MEG'). 
%          'decoding_window' (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%           other possible name-value pairs: SVM settings, svm evaluation metrics (see Documentation).
% Output: time and space-resolved accuracy and F1-score (channels x time).
% Performs time and space-resolved SVM decoding of MEG data, using stratified k-fold cross-validation and LibLinear for each time window and sensor cluster.
% Uses sensor neighbours defined using FT template configuration (or custom structure containing source indices).

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

if size(train_data,1)~=size(test_data,1) || size(train_data,2)~=size(test_data,2)
    error('Training and test data need to have equal numbers of features')
end;

parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','decoding_window', 'time', 'pseudo'}); %converted struct will be fed into decoding function
clear p;

%create time axis
if ~isempty(dec_args.time)
    time = dec_args.time;
else
    time = 1:size(train_data,2);
end;

if length(time)~=size(train_data,2)
    time = 1:size(train_data,2);
    fprintf('Warning: time axis does not match dataset size. Replacing with default time axis...');
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
        lims(2) = size(train_data,2);
    end;
        
else
    lims = [1 size(train_data,2)];
end;

if size(train_data,2)<lims(2)
    lims(2) = size(train_data,2);
end;


fprintf('\nRunning classifier... '); 
%loop through time
train_data_svm = arrayfun(@(i) reshape(train_data(:, i:i+dec_args.window_length-1,:), size(train_data,1)*dec_args.window_length, size(train_data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %time selection
test_data_svm = arrayfun(@(i) reshape(test_data(:, i:i+dec_args.window_length-1,:), size(test_data,1)*dec_args.window_length, size(test_data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %time selection
results_tmp = arrayfun(@(i) svm_decode_holdout(train_data_svm{i},train_labels, test_data_svm{i}, test_labels, svm_par), 1:length(train_data_svm));
results.Accuracy = cell2mat({results_tmp.Accuracy});
results.WeightedFscore = cell2mat({results_tmp.WeightedFscore});
if svm_par.weights
    results.Weights =  cat(1,results_tmp(:).Weights)';
    results.WeightsPatterns =  cell2mat({results_tmp.WeightPatterns});
end
results.Confusion = cat(3,results_tmp(:).Confusion);
results.Sensitivity = cell2mat({results_tmp(:).Sensitivity});
results.Specificity = cell2mat({results_tmp(:).Specificity});
clear train_data_svm test_data_svm;

end

