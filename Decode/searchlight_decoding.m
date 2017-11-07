function [ accuracy, results] = searchlight_decoding( data, labels, info_file, varargin )
% Inputs: data, labels, info_file (.mat file obtained using get_sensor_info or get_source_info).
% Optional: channel set (string or  cell array of strings; default: 'MEG'), 
% decoding window (limits; default: [-0.1 0.9])
% and window length (in sampled time points; default here:10);
% SVM settings, svm evaluation metrics (see Documentation).
% Output: time and space-resolved accuracy.
% Performs time and space-resolved SVM decoding of MEG data. Uses
% sensor neighbours defined using FT template configuration (or custom structure).

%parse inputs
dec_args = decoding_args;
svm_par = svm_args;
svm_results = svm_eval;
list = [fieldnames(dec_args); fieldnames(svm_par); fieldnames(svm_results)];
p = inputParser;

for i = 1:length(properties(decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
for ii = i+1:length(properties(dec_args))+length(properties(svm_args))
    addParameter(p, list{ii}, svm_par.(list{ii}));
end;
for i = ii+1:length(properties(dec_args))+length(properties(svm_args))+length(properties(svm_eval))
    addParameter(p, list{i}, svm_results.(list{i}));
end;

parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','decoding_window'}); %converted struct will be fed into decoding function
clear p svm_results;

load (info_file);

%time limits for decoding window
if ~isempty(dec_args.decoding_window)
    lims = [find(round(time,3)==dec_args.decoding_window(1)) find(round(time,3)==dec_args.decoding_window(2))];
else
    lims = [1 length(time)];
end;

%channel indices for each searchlight iteration
if exist ('chan_labels', 'var') %the sensor-space case
    chan_idx = arrayfun(@(i) find(ismember(chan_labels,neighbours(i).neighblabel)), 1:length(neighbours), 'UniformOutput', false); %store all searchlight idx in a cell array
else
    chan_idx = source_idx; %the source-space case
end;
accuracy = zeros(length(chan_idx), length(lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1));
fprintf('\nRunning searchlight '); 

results = cell(1,length(chan_idx));

for c = 1:length(chan_idx)
    
    fprintf('%d out of %d', c, length(chan_idx));
    data_svm = arrayfun(@(i) reshape(data(chan_idx{c}, i:i+dec_args.window_length-1,:), length(chan_idx{c})*dec_args.window_length, size(data,3))', lims(1):dec_args.window_length:lims(2)-dec_args.window_length+1, 'UniformOutput', false); %channel and time selection
    results = arrayfun(@(i) svm_decode_kfold(data_svm{i},labels,svm_par), 1:length(data_svm));
    accuracy(c,:) = cell2mat({results.Accuracy});
    fprintf((repmat('\b',1,numel([num2str(c-1) num2str(length(chan_idx))])+8)));
    
end;

end

