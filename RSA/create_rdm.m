function [ rdm, rdm_tril ] = create_rdm( data, varargin)
%Creates MEG RDM from trial x channel x time data.
%Input: data + optional arguments in name-value pairs (channel, time window, window length)
%Output: Time-resolved MEG RDM in matrix and vector formats.

opt = decoding_args;
list = fieldnames(opt);
p = inputParser;
for i = 1:length(list)
    addParameter(p, list{i}, opt.(list{i}));
end;
addParameter(p, 'sensor_idx', []);
addParameter(p, 'normalize', 2); %0: don't normalize; 1:normalize per time point; 2:across time points
parse(p, varargin{:});
args = p.Results;
clear p opt;

%create time axis
if ~isempty(args.time)
    time = args.time;
elseif ~isempty(args.sensor_idx) && isfield(args.sensor_idx, 'time')
    time = args.sensor_idx.time;
else
    time = 1:size(data,ndims(data));
end;


%time limits for decoding window
if ~isempty(args.decoding_window)
    lims = [find(round(time,3)==args.decoding_window(1)) find(round(time,3)==args.decoding_window(2))];
else
    lims = [1 length(time)];
end;

%here can simply include all possibilities of spatial selection 
%i.e. sensor sets, searchlight, AAL & source indices
if strcmp(args.channels, 'searchlight')
    if isstruct (args.sensor_idx) %the sensor-space case
        chan_idx = arrayfun(@(i) find(ismember({sensor_idx.label},[sensor_idx(i).label; sensor_idx(i).neighblabel])), 1:length(sensor_idx), 'UniformOutput', false); %store all searchlight idx in a cell array
    else
        chan_idx = source_idx; %the source-space case
    end;
    
    %here write the searchlight itself
    tp = lims(1):args.window_length:lims(2)-args.window_length+1;
    
    rdm = zeros(size(data,1), size(data,1), length(chan_idx), length(tp));
    rdm_tril = zeros((size(data,1)^2)/2 - size(data,1)/2, length(chan_idx), length(tp)); %lower triangular part of matrix for correlations
    
    for c = 1:length(chan_idx)
        
        for t = 1:length(tp)
            D = squareform(pdist(reshape(data(:,chan_idx{c},tp(t):tp(t)+args.window_length-1), size(data,1), length(chan_idx{c})*args.window_length)));
            D = (D-min(D(:)))./(max(D(:))-min(D(:)));
            rdm(:,:,c,t) = D;
            rdm_tril(:,c,t) = D(tril(true(size(D)),-1));
        end;
    
    end;
    
    
else
    chan_idx = 1:size(data,2);
    if ~strcmp (args.channels, 'MEG') %if we need to subselect sensors
        chan = [];
        for i = 1:length(args.channels)
            idx = cellfun('isempty',strfind({sensor_idx.label},args.channels{i}));
                chan = [chan chan_idx(~idx)]; %#ok<AGROW>
        end;
        chan_idx = chan;
    end;
    
    data = data(:,chan_idx, lims(1):lims(2));
    tp = lims(1):args.window_length:lims(2)-args.window_length+1;
    
    rdm = zeros(size(data,1), size(data,1),length(tp));
    rdm_tril = zeros((size(data,1)^2)/2 - size(data,1)/2, length(tp));
    
    for t = 1:length(tp)
        D = squareform(pdist(reshape(data(:,:,tp(t):tp(t)+args.window_length-1), size(data,1), size(data,2)*args.window_length)));
        if args.normalize==1
            D = (D-min(D(:)))./(max(D(:))-min(D(:)));
        end;
        rdm(:,:,t) = D;
        rdm_tril(:,t) = D(tril(true(size(D)),-1));
    end;
    
    %normalize across time points
    if args.normalize==2
        rdm = (rdm-min(rdm(:)))./(max(rdm(:))-min(rdm(:)));
        rdm_tril = (rdm_tril-min(rdm_tril(:)))./(max(rdm_tril(:))-min(rdm_tril(:)));
    end;
    
end;

    
end

