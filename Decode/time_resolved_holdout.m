function [ results ] = time_resolved_holdout( train_data, train_labels, test_data, test_labels, varargin )
% Inputs: training data, training labels, test data, test labels.
% Data format: channels/sources x time x trials.
% Optional: 
%          'channels', channel set (string or cell array of strings; default: 'MEG').
%          'decoding_window' (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%          'pseudo', default [], create pseudotrials: example [5 100], average groups of 5 trials with 100 random assignments of trials to groups
%          'mnn', default true, perform multivariate noise normalization (recommended, Guggenmos et al. 2018)
%
%           Below are other name-value pairs that control the SVM settings, with defaults:
%           
%          solver = 1; %only applies to liblinear: 1: L2 dual-problem; 2: L2 primal; 3:L2RL1L...
%          boxconstraint = 1; --> C-parameter: Note, we don't have any options for optimizing this, need to write it separately if needed
%          standardize = true; --> standardize features using mean and SD of training set (recommended)
%          weights = false; --> calculate weights (by retraining model on whole dataset)
%
% Outputs: structure containing classification performance metrics (for each timepoint and cross-validation round).
%
% DC Dima 2018 (diana.c.dima@gmail.com)

%parse inputs
dec_args = args.decoding_args;
svm_par = args.svm_args;
list = [fieldnames(dec_args); fieldnames(svm_par)];
p = inputParser;

for i = 1:length(properties(args.decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
for ii = i+1:length(properties(args.decoding_args))+length(properties(args.svm_args))
    addParameter(p, list{ii}, svm_par.(list{ii}));
end;

if size(train_data,1)~=size(test_data,1) || size(train_data,2)~=size(test_data,2)
    error('Training and test data need to have equal numbers of features')
end;

parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','decoding_window', 'time', 'pseudo','mnn'}); %converted struct will be fed into decoding function
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

%create pseudo-trials if requested
if ~isempty(dec_args.pseudo)
    data_tmp = cat(3,train_data,test_data);
    labels_tmp = [train_labels;test_labels];
    [data,labels] = create_pseudotrials(data_tmp, labels_tmp, dec_args.pseudo(1), dec_args.pseudo(2));
    data = cat(3,data{:}); labels = cat(1,labels{:});
    train_data = data(:,:,1:length(train_labels)); 
    train_labels = labels(1:length(train_labels));
    test_data = data(:,:,(length(train_labels)+1):(length(train_labels)+length(test_labels)));
    test_labels = labels(length(train_labels)+1):(length(train_labels)+length(test_labels));
    clear data labels data_tmp labels_tmp;
end

%whiten data if requested
if dec_args.mnn
    [train_data,test_data] = whiten_data(train_data,train_labels,test_data);
end

fprintf('\nRunning classifier... '); 
%loop through time
train_data_svm = arrayfun(@(i) reshape(train_data(:, i:i+dec_args.window_length-1,:), size(train_data,1)*dec_args.window_length, size(train_data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %time selection
test_data_svm = arrayfun(@(i) reshape(test_data(:, i:i+dec_args.window_length-1,:), size(test_data,1)*dec_args.window_length, size(test_data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %time selection
results_tmp = arrayfun(@(i) svm_decode_holdout(train_data_svm{i},train_labels, test_data_svm{i}, test_labels, svm_par), 1:length(train_data_svm));
results.Accuracy = cell2mat({results_tmp.Accuracy});
results.WeightedFscore = cell2mat({results_tmp.WeightedFscore});
if svm_par.weights
    results.Weights =  cat(1,results_tmp(:).Weights)';
    results.WeightPatterns =  cell2mat({results_tmp.WeightPatterns});
    results.WeightPatternsNorm =  cell2mat({results_tmp.WeightPatternsNorm});
    if dec_args.window_length>1
        results.Weights = reshape(results.Weights, size(train_data,1), dec_args.window_length, size(results.Weights,2));
        results.Weights = squeeze(mean(results.Weights,2));
        results.WeightPatterns = reshape(results.WeightPatterns, size(train_data,1), dec_args.window_length, size(results.WeightPatterns,2));
        results.WeightPatterns = squeeze(mean(results.WeightPatterns,2));
end
results.Confusion = cat(3,results_tmp(:).Confusion);
results.Sensitivity = cell2mat({results_tmp(:).Sensitivity});
results.Specificity = cell2mat({results_tmp(:).Specificity});
results.PredictedLabels = cell2mat({results_tmp(:).PredictedLabels});
if ~isempty(dec_args.pseudo)
    results.PseudoLabels = test_labels;
end
clear train_data_svm test_data_svm;

end

