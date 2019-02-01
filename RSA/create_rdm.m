function [ rdm, rdm_tril ] = create_rdm( data, varargin)
% Creates Euclidean distance MEG RDM from trial x channel x time data. Can create a time-resolved or spatiotemporally resolved RDM
% Input: data (trial x channel/source x time)
% Optional arguments in name-value pairs: 
%          'channels', channel set (string or cell array of strings; default: 'MEG').
%          'decoding_window' in this case, RSA analysis window (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%          'pseudo', default [], create pseudotrials: example [5 100], average groups of 5 trials with 100 random assignments of trials to groups
%          'sensor_idx', default [], sensor neighbours structure or cell array with source space searchlights/ROIs
%          'normalize', default 0, min-max RDM normalization: 0:don't normalize; 1:normalize per time point; 2:normalize across time points
%
%Output: Time-resolved MEG RDM in matrix and vector formats. (Feature x space x time, OR feature x time)
%
%DC Dima 2017 (diana.c.dima@gmail.com)

addpath(strrep(mfilename('fullpath'), fullfile('RSA','create_rdm'),'Decode'));

opt = args.decoding_args;
list = fieldnames(opt);
p = inputParser;
for i = 1:length(list)
    addParameter(p, list{i}, opt.(list{i}));
end;
addParameter(p, 'sensor_idx', []);
addParameter(p, 'normalize', 0); 
parse(p, varargin{:});
rsa_args = p.Results;
clear p opt;

%create time axis
if ~isempty(rsa_args.time)
    time = rsa_args.time;
elseif ~isempty(rsa_args.sensor_idx) && isfield(rsa_args.sensor_idx, 'time')
    time = rsa_args.sensor_idx.time;
else
    time = 1:size(data,ndims(data));
end;


%time limits for analysis window
if ~isempty(rsa_args.decoding_window)
    lims = [find(round(time,3)==rsa_args.decoding_window(1)) find(round(time,3)==rsa_args.decoding_window(2))];
else
    lims = [1 length(time)];
end;

%here can simply include all possibilities of spatial selection 
%i.e. sensor sets, searchlight, AAL & source indices
if ~isempty(rsa_args.sensor_idx)
    if isstruct (rsa_args.sensor_idx) %the sensor-space case
        chan_idx = arrayfun(@(i) find(ismember({sensor_idx.label},[sensor_idx(i).label; sensor_idx(i).neighblabel])), 1:length(sensor_idx), 'UniformOutput', false); %store all searchlight idx in a cell array
    else
        chan_idx = source_idx; %the source-space case
    end;
    
    %here write the searchlight itself
    tp = lims(1):rsa_args.window_length:lims(2)-rsa_args.window_length+1;
    
    rdm = zeros(size(data,1), size(data,1), length(chan_idx), length(tp));
    rdm_tril = zeros((size(data,1)^2)/2 - size(data,1)/2, length(chan_idx), length(tp)); %lower triangular part of matrix for correlations
    
    for c = 1:length(chan_idx)
        
        for t = 1:length(tp)
            D = squareform(pdist(reshape(data(:,chan_idx{c},tp(t):tp(t)+rsa_args.window_length-1), size(data,1), length(chan_idx{c})*rsa_args.window_length)));
            D = (D-min(D(:)))./(max(D(:))-min(D(:)));
            rdm(:,:,c,t) = D;
            rdm_tril(:,c,t) = D(tril(true(size(D)),-1));
        end;
    
    end;
    
    
else
    chan_idx = 1:size(data,2);
    if ~strcmp (rsa_args.channels, 'MEG') %if we need to subselect sensors
        chan = [];
        for i = 1:length(rsa_args.channels)
            idx = cellfun('isempty',strfind({sensor_idx.label},rsa_args.channels{i}));
                chan = [chan chan_idx(~idx)]; %#ok<AGROW>
        end;
        chan_idx = chan;
    end;
    
    data = data(:,chan_idx, lims(1):lims(2));
    tp = lims(1):rsa_args.window_length:lims(2)-rsa_args.window_length+1;
    
    rdm = zeros(size(data,1), size(data,1),length(tp));
    rdm_tril = zeros((size(data,1)^2)/2 - size(data,1)/2, length(tp));
    
    for t = 1:length(tp)
        D = squareform(pdist(reshape(data(:,:,tp(t):tp(t)+rsa_args.window_length-1), size(data,1), size(data,2)*rsa_args.window_length)));
        if rsa_args.normalize==1
            D = (D-min(D(:)))./(max(D(:))-min(D(:)));
        end;
        rdm(:,:,t) = D;
        rdm_tril(:,t) = D(tril(true(size(D)),-1));
    end;
    
    %normalize across time points
    if rsa_args.normalize==2
        rdm = (rdm-min(rdm(:)))./(max(rdm(:))-min(rdm(:)));
        rdm_tril = (rdm_tril-min(rdm_tril(:)))./(max(rdm_tril(:))-min(rdm_tril(:)));
    end;
    
end;

    
end

