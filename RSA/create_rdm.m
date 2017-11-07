function [ rdm, rdm_triu ] = create_rdm( data, neighbours, varargin)
%Creates MEG RDM from channel x time x trial data.
%Input: data + optional arguments in name-value pairs (channel, time window, window length)
%Output: Time-resolved MEG RDM in matrix and vector formats.

opt = lcmv_args;
list = fieldnames(opt);
p = inputParser;
for i = 1:length(list)
    addParameter(p, list{i}, opt.(list{i}));
end;
parse(p, varargin{:});
dec_args = p.Results;
clear p opt;

load (neighbours);

%time limits for decoding window
if ~isempty(dec_args.decoding_window)
    lims = [find(round(time,3)==dec_args.decoding_window(1)) find(round(time,3)==dec_args.decoding_window(2))];
else
    lims = [1 length(time)];
end;

%here can simply include all possibilities of spatial selection 
%i.e. sensor sets, searchlight, AAL & source indices
if strcmp(dec_args.channels, 'searchlight')
    if exist ('neighbours', 'var') %the sensor-space case
        chan_idx = arrayfun(@(i) find(ismember(chan_labels,neighbours(i).neighblabel)), 1:length(neighbours), 'UniformOutput', false); %store all searchlight idx in a cell array
    else
        chan_idx = source_idx; %the source-space case
    end;
    
    %here write the searchlight itself
    tp = lims(1):opt.window_length:lims(2);
    
    rdm = zeros(size(data,3), size(data,3),length(tp), length(chan_idx));
    rdm_triu = zeros((size(data,3)^2)/2 + size(data,3)/2, length(tp), length(chan_idx));
    
    for c = 1:length(chan_idx)
        
        for t = 1:length(tp)
            D = squareform(pdist(reshape(data(chan_idx{c},tp(t):tp(t)+opt.window_length-1,:), size(data,3), length(chan_idx{c})*opt.window_length)));
            D = (D-min(D(:)))./(max(D(:))-min(D(:)));
            rdm(:,:,t,c) = D;
            rdm_triu(:,t,c) = D(triu(true(size(D))));
        end;
    
    end;
    
    
else
    chan_idx = 1:length(chan_labels);
    if ~strcmp (dec_args.channels, 'MEG') %if we need to subselect sensors
        chan = [];
        for i = 1:length(dec_args.channels)
            idx = cellfun('isempty',strfind(chan_labels,dec_args.channels{i}));
            chan = [chan chan_idx(~idx)]; %#ok<AGROW>
        end;
        chan_idx = chan;
    end;
    
    data = data(chan_idx, lims(1):lims(2), :);
    tp = lims(1):opt.window_length:lims(2);
    
    rdm = zeros(size(data,3), size(data,3),length(tp));
    rdm_triu = zeros((size(data,3)^2)/2 + size(data,3)/2, length(tp));
    
    for t = 1:length(tp)
        D = squareform(pdist(reshape(data(:,tp(t):tp(t)+opt.window_length-1,:), size(data,3), size(data,1)*opt.window_length)));
        D = (D-min(D(:)))./(max(D(:))-min(D(:)));
        rdm(:,:,t) = D;
        rdm_triu(:,t) = D(triu(true(size(D))));
    end;
    
end;

    
end

