function [data] = read_meg_data(dataset, condition, varargin)
% Inputs: dataset, condition; optional: prestimulus, poststimulus,
% baseline, bandpass filter, resampling frequency
% Output: preprocessed MEG data. Channels x time x trials matrix.
% Uses fieldtrip toolbox (Oostenveld et al., 2011) to read in MEG data
% for a specified condition.

opt = preproc_args;
list = fieldnames(opt);
p = inputParser;
for i = 1:length(list)
    addParameter(p, list{i}, opt.(list{i}));
end;
parse(p, varargin{:});

cfg = [];
cfg.dataset = dataset;

%read in data
cfg.trialdef.prestim = opt.prestim;
cfg.trialdef.poststim = opt.poststim;
cfg.trialdef.eventtype = condition;
cfg.trialfun = 'ft_trialfun_NoBadTrials';
cfg = ft_definetrial(cfg);

%preprocessing options
cfg.dftfilter = 'yes';            %remove mains noise
if ~isempty(opt.bandpass_filter)
    if opt.bandpass_filter(1) == 0
        cfg.lpfilter = 'yes';
        cfg.lpfreq = opt.bandpass_filter(2);
    elseif opt.bandpass_filter(2) == 0
        cfg.hpfilter = 'yes';
        cfg.hpreq = opt.bandpass_filter(1);
        cfg.padding = prestim + poststim + 4; %2 s of padding both sides
        cfg.padtype = 'mirror';
    else
        cfg.bpfilter = 'yes';
        cfg.bpfreq = opt.bandpass_filter;
        cfg.bpfiltord = 3;
        cfg.padding = prestim + poststim + 4; %2 s of padding both sides
        cfg.padtype = 'mirror';
    end;
end;

cfg.channel = 'MEG';

if ~isempty(opt.baseline)
    cfg.demean = 'yes';
    cfg.baselinewindow = opt.baseline; %baselining
    cfg.detrend = 'yes';
end;
data = ft_preprocessing(cfg);

%here we recover original time axis after padding
if isfield(cfg, 'padding')
    for i = 1:length(data.time)
        data.time{i} = data.time{i}-2;
    end;
end;

if ~isempty(opt.resamplefs)
    cfg = [];
    cfg.resamplefs = opt.resamplefs;
    cfg.detrend='no';
    data=ft_resampledata(cfg,data);
end;

data = cat(3,data.trial{:});

end

