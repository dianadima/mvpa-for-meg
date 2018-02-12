function [ results ] = temp_gen_decoding( data, labels, varargin )
% Performs temporal generalization SVM decoding of MEG data (train on each time point/window, test on all others).
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
if ~isempty(dec_args.sensor_idx)
    if ~exist('chan_idx', 'var')
        chan_idx = 1:size(data,1); %initialize with entire array
        if ~strcmp (dec_args.channels, 'MEG') %if we need to subselect sensors
            chan = [];
            for i = 1:length(dec_args.channels)
                idx = cellfun('isempty',strfind({sensor_idx.label},dec_args.channels{i}));
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
data_svm = arrayfun(@(i) reshape(data(chan_idx,i:i+dec_args.window_length-1,:), length(chan_idx)*dec_args.window_length, size(data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %features are channelsxtime,observations are trials    
results = zeros(length(data_svm), length(data_svm));

%and here we run the classifier -first train on all time points and store the models
svm_model = arrayfun(@(t) svm_train(data_svm{t}, labels, svm_par), 1:length(data_svm), 'UniformOutput', false);
fprintf('\nFinished training all models...\r')  
fprintf('\rNow testing...\r')

%testing on all time points
for t = 1:length(svm_model)
    
    if mod(t,50)==0
        fprintf('\n%d out of %d...\n', t, length(svm_model))
    end;
    
    for i = 1:length(data_svm)
        
        %standardize test data, using values based on training data
        if svm_par.standardize
            test_data = (data_svm{i} - repmat(min(data_svm{t}, [], 1), size(data_svm{i}, 1), 1)) ./ repmat(max(data_svm{t}, [], 1) - min(data_svm{t}, [], 1), size(data_svm{i}, 1), 1);
        end;
        
        [~, accuracy, ~] = predict(labels, sparse(test_data), svm_model{t}, '-q 1');
        results(t,i) = accuracy(1);

    end;
    
end;
            
end

