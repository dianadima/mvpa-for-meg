function [ ] = plot_searchlight_results(results, neighbours, varargin)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

dec_args = decoding_args;
dec_args.window_length = 0.05; %just this one different than default
list = fieldnames(dec_args);
p = inputParser;
for i = 1:length(properties(decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
addParameter(p, 'clim', [40 100]);
addParameter(p, 'colormap', 'jet');
addParameter(p, 'highlight_channels', []);
addParameter(p, 'configuration', 'CTF275');
parse(p, varargin{:});
dec_args = p.Results;

%create time axis
if ~isempty(dec_args.time)
    time = dec_args.time;
elseif isfield(neighbours, 'time')
    time = neighbours.time;
else
    if ismatrix(results)
        time = 1:size(results,2);
    elseif ndims(results)==3
        time = 1:size(results,3);
    end;
end;

if ~isempty(dec_args.decoding_window)
    if ~isempty(find(round(time,3)==dec_args.decoding_window(1),1))
        lims(1) = find(round(time,3)==dec_args.decoding_window(1));
    else
        fprintf('Warning: starting timepoint not found. Starting from 1...');
        lims(1) = 1;
    end;
    if ~isempty(find(round(time,3)==dec_args.decoding_window(2),1))
        lims(2) = find(round(time,3)==dec_args.decoding_window(2));
    else
        fprintf('Warning: end timepoint not found. Plotting til the end...');
        lims(2) = length(time);
    end;
        
else
    lims = [1 length(time)];
end;

if length(time)<lims(2)
    lims(2) = length(time);
end;

if ismatrix(results)
    acc = results;
elseif ndims(results)==3
    acc = squeeze(mean(results,1));
    fprintf('Warning: assuming subjects are 1st dimension of accuracy matrix....')
else
    error('Results should be a 2d or 3d matrix containing subjects x channels x time');
end;

%create structure for plotting
sens.acc = acc;
sens.label = {neighbours(:).label};
sens.dimord = 'chan_time';
sens.time = time;
cfg.xlim = time(lims(1)):dec_args.window_length:time(lims(2));
%1D data case (1 timepoint)
if numel(cfg.xlim)==1
    cfg.xlim(2) = cfg.xlim(1);
end;    
cfg.zlim = dec_args.clim;
cfg.layout = [dec_args.configuration '.lay'];
cfg.parameter = 'acc';
cfg.style = 'straight';
cfg.comment = 'no';
cfg.interactive = 'no';
cfg.colormap = dec_args.colormap;
if ~isempty(dec_args.highlight_channels)
    cfg.highlight = 'on';
    cfg.highlightchannel = dec_args.highlight_channels;
    cfg.highlightmarker = 'o';
    cfg.highlightcolor = [1 1 1];
    cfg.highlightsize = 8;
end;

ft_topoplotER(cfg,sens)

end

