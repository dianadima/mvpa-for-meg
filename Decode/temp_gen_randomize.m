function [ results ] = temp_gen_randomize( data, labels, varargin )
% Performs temporal generalization SVM decoding of MEG data (train on each time point/window, test on all others). Uses hold-out validation (train on half, test on half).
% Inputs: data, labels.
% Optional: 'sensor_idx', structure obtained using get_sensor_info - for channel selection. You can also just provide numerical indices, in which case you don't need the structure. 
%                        If you want to subselect features on source space data, you need to manually provide numerical indice (i.e. 'channels', [1:1000]). 
%          'channels', channel set set (string or cell array of strings; default: 'MEG'). 
%          'decoding_window' (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%                   other possible name-value pairs: SVM settings, svm evaluation metrics (see Documentation).
%
% Outputs: time x time accuracy matrix.


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
addParameter(p, 'sensor_idx', []);
parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','decoding_window', 'time', 'sensor_idx'}); %converted struct will be fed into decoding function
clear p;

%get channel indices and time axis. Numerical channel indices take priority
if ~iscell(dec_args.channels) && ~ischar(dec_args.channels)
    chan_idx = dec_args.channels;
end;
if (~isempty(dec_args.sensor_idx)) && (~strcmp (dec_args.channels, 'MEG'))
    if ~exist('chan_idx', 'var')
        chan = [];
        for i = 1:length(dec_args.channels)
            idx = cellfun('isempty',strfind({sensor_idx.label},dec_args.channels{i}));
            chan = [chan chan_idx(~idx)]; %#ok<AGROW>
        end;
        chan_idx = chan;
    end;
else
    if ~exist('chan_idx', 'var')
        chan_idx = 1:size(data,1);
    end;
end;

%create time axis
if ~isempty(dec_args.time)
    time = dec_args.time;
elseif ~isempty(dec_args.sensor_idx) && isfield(sensor_idx, 'time')
    time = sensor_idx.time;
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
%we use a hold-out method that trains on half and tests on the other half, for speed reasons; other methods such as nested kfold can be used
classes = unique(labels); idx1 = find(labels==classes(1)); idx2 = find(labels==classes(2));
train_idx = [idx1(randperm(floor(length(idx1)/2))); idx2(randperm(floor(length(idx2)/2)))];
test_idx = [idx1; idx2]; test_idx(train_idx) = [];

train_data = arrayfun(@(i) reshape(data(chan_idx,i:i+dec_args.window_length-1,train_idx), length(chan_idx)*dec_args.window_length, length(train_idx))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %features are channelsxtime,observations are trials    
test_data = arrayfun(@(i) reshape(data(chan_idx,i:i+dec_args.window_length-1,test_idx), length(chan_idx)*dec_args.window_length, length(test_idx))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %features are channelsxtime,observations are trials    

results = zeros(length(train_data), length(train_data));

%and here we run the classifier -first train on all time points and store the models
labelst = labels(train_idx); labelst = labelst(randperm(length(labelst)));
svm_model = arrayfun(@(t) svm_train(train_data{t}, labelst, svm_par), 1:length(train_data), 'UniformOutput', false);
fprintf('\nFinished training all models...\r')  
fprintf('\rNow testing...\r')

%testing on all time points
for t = 1:length(svm_model)
    
    if mod(t,50)==0
        fprintf('\n%d out of %d...\n', t, length(svm_model))
    end;
    
    labelst = labels(test_idx); labelst = labelst(randperm(length(labelst)));
    
    for i = 1:length(test_data)
        
        %standardize test data, using values based on training data
        if svm_par.standardize
            test_data_i = (test_data{i} - repmat(min(train_data{t}, [], 1), size(test_data{i}, 1), 1)) ./ repmat(max(train_data{t}, [], 1) - min(train_data{t}, [], 1), size(test_data{i}, 1), 1);
        end;    
        
        [~, accuracy, ~] = predict(labelst, sparse(test_data_i), svm_model{t}, '-q 1');
        results(t,i) = accuracy(1);

    end;
    
end;
            
end
