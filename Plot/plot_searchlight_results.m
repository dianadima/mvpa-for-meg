function [ ] = plot_searchlight_results(results, info_file, varargin)
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
parse(p, varargin{:});
dec_args = p.Results;

load(info_file);

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
    error('Results should be a 2d or 3d matrix containing subjects x time x channels');
end;

%create structure for plotting
sens.acc = acc';
sens.label = chan_labels;
sens.dimord = 'chan_time';
sens.time = time(lims(1)):dec_args.window_length:time(lims(2));
cfg.xlim = time(lims(1)):dec_args.window_length:time(lims(2));
cfg.zlim = dec_args.clim;
cfg.layout = 'CTF275.lay';
cfg.parameter = 'acc';
cfg.style = 'straight';
ft_topoplotER(cfg,sens)

end

