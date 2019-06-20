function [ accuracy, Fscore] = searchlight_holdout( train_data, train_labels, test_data, test_labels, cluster_idx, varargin )
% Inputs: training data, training labels, test data, test labels. Data format: channels/sources x time x trials.
%         cluster_idx, structure obtained using get_sensor_info or get_source_info - for channel/source selection. You can also just provide numerical indices, in which case you don't need the structure.
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
% Outputs: accuracy, F-score (clusters x time)
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
clear p

%channel indices for each searchlight iteration
if isstruct(cluster_idx) %the sensor-space case
    chan_idx = arrayfun(@(i) find(ismember({cluster_idx.label},[cluster_idx(i).label; cluster_idx(i).neighblabel])), 1:length(cluster_idx), 'UniformOutput', false); %store all searchlight idx in a cell array
else
    chan_idx = cluster_idx; %the source-space case
end

%create time axis
if ~isempty(dec_args.time)
    time = dec_args.time;
elseif isfield(cluster_idx, 'time')
    time = cluster_idx(1).time;
else
    time = 1:size(train_data,2);
end

if length(time)~=size(train_data,2)
    time = 1:size(train_data,2);
    fprintf('\nWarning: time axis does not match dataset size. Replacing with default time axis...');
end

%time limits for decoding window
if ~isempty(dec_args.decoding_window)
    lims(1) = nearest(time,dec_args.decoding_window(1));
    lims(2) = nearest(time,dec_args.decoding_window(2));       
else
    lims = [1 size(train_data,2)];
end

if size(train_data,2)<lims(2)
    lims(2) = size(train_data,2);
end;

%create pseudo-trials if requested - separately for training & test data
if ~isempty(dec_args.pseudo)
    [train_data,train_labels] = create_pseudotrials(train_data, train_labels, dec_args.pseudo(1), dec_args.pseudo(2));
    [test_data,test_labels] = create_pseudotrials(test_data, test_labels, dec_args.pseudo(1), dec_args.pseudo(2));
end

%whiten data if requested
if dec_args.mnn
    [train_data,test_data] = whiten_data(train_data,train_labels,test_data);
end

accuracy = zeros(length(chan_idx), length(lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1));
Fscore = zeros(length(chan_idx), length(lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1));
fprintf('\nRunning searchlight '); 

%loop through channels and time
for c = 1:length(chan_idx)
    
    fprintf('%d out of %d', c, length(chan_idx));
    train_data_svm = arrayfun(@(i) reshape(train_data(chan_idx{c}, i:i+dec_args.window_length-1,:), length(chan_idx{c})*dec_args.window_length, size(train_data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %channel and time selection
    test_data_svm = arrayfun(@(i) reshape(test_data(chan_idx{c}, i:i+dec_args.window_length-1,:), length(chan_idx{c})*dec_args.window_length, size(test_data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %channel and time selection
    results = arrayfun(@(i) svm_decode_holdout(train_data_svm{i},train_labels, test_data_svm{i}, test_labels, svm_par), 1:length(train_data_svm));
    accuracy(c,:) = cell2mat({results.Accuracy});
    Fscore(c,:) = cell2mat({results.WeightedFscore});
    fprintf((repmat('\b',1,numel([num2str(c) num2str(length(chan_idx))])+8)));
    clear train_data_svm test_data_svm;
    
end;

end

