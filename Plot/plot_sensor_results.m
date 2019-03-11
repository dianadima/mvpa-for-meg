function [ ] = plot_sensor_results(results, neighbours, varargin)
% Plots results of sensor-space searchlight decoding.
% Inputs: results: matrix of accuracy/decoding performance. Must be channels, channels x time, or subjects x channels x time.
%         Note that if you specify a 'window_length'>1, accuracy will be averaged over the time window (because it uses ft_topoplotER). 
%         neighbours: sensor grouping structure obtained using get_sensor_info (i.e., fieldtrip function prepare_neighbours).
% Optional inputs:
%   'configuration' (default 'CTF275'): specify configuration of sensors. The fieldtrip layout <configuration>.lay will be loaded for plotting.
%   'colorlim' (default [40 100]): colour limits
%   'colormap' (default 'jet')
%   'highlight_channels', significant channels can be marked with white asterisks. Cell array of strings (e.g., {'MRO18', 'MRO19'}.
%   'time', time axis, if you wish to specify time units, rather than sampled points.
%   'decoding_window', which portion of the data is being plotted.
%   'window_length', what size time window is included in one plot. This may lead to multiple plots. With a decoding window of [0 1] and a
%                   window length of 0.5, you will get two plots showing 50 ms each (time axis will become 0:0.5:1).
%
% DC Dima 2018 (diana.c.dima@gmail.com)

dec_args = decoding_args;
dec_args.window_length = 0.05; %just this one different than default
list = fieldnames(dec_args);
p = inputParser;
for i = 1:length(properties(decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
addParameter(p, 'colorlim', [40 100]);
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
    if ~isempty(find(round(time,3)==dec_args.decoding_window(end),1))
        lims(2) = find(round(time,3)==dec_args.decoding_window(end));
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
cfg.zlim = dec_args.colorlim;
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

