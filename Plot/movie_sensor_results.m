function [ ] = movie_sensor_results( results, neighbours, output_file, varargin )
% Plot sensor-space searchlight decoding results as a movie.
% Inputs: results: matrix of accuracy/decoding performance. Must be channels x time, or subjects x channels x time.
%         neighbours: sensor grouping structure obtained using get_sensor_info (i.e., fieldtrip function prepare_neighbours).
%         output_file: movie filename
% Optional inputs:
%   'configuration' (default 'CTF275'): specify configuration of sensors. The fieldtrip layout <configuration>.lay will be loaded for plotting.
%   'colorlim' (default [40 100]): colour limits
%   'colormap' (default 'jet')
%   'result_type' (default 'Accuracy (%)'): will be plotted as colorbar axis
%   'visible' (default 'on'): show figure
%
% DC Dima 2018 (diana.c.dima@gmail.com)

if isempty(neighbours)
    [~, ftdir] = ft_version; %get FT directory
    load([ftdir '/template/neighbours/ctf275_neighb.mat']);
end;

dec_args = args.decoding_args;
dec_args.decoding_window = []; %replace default
list = fieldnames(dec_args);
p = inputParser;
for i = 1:length(properties(args.decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
addParameter(p, 'configuration', 'CTF275');
addParameter(p, 'colorlim', [40 100]);
addParameter(p, 'colormap', 'jet');
addParameter(p, 'result_type', 'Accuracy (%)');
addParameter(p, 'visible', 'on');
parse(p, varargin{:});
dec_args = p.Results;

%create time axis
if ~isempty(dec_args.time)
    time = dec_args.time;
elseif isfield(neighbours, 'time') && ismember(length(neighbours(1).time), [size(results,2), size(results,3)])
    time = neighbours(1).time;
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
acc = acc(:,lims(1):lims(2));
sens.label = {neighbours(:).label};
sens.dimord = 'chan_time';
sens.time = 1;
cfg.xlim = [1 1];
cfg.zlim = dec_args.colorlim;
cfg.layout = [dec_args.configuration '.lay'];
cfg.parameter = 'acc';
cfg.style = 'straight';
cfg.comment = 'no';
cfg.interactive = 'no';
cfg.colormap = dec_args.colormap;
cfg.visible = p.Results.visible;

F(size(acc,2)) = struct('cdata',[],'colormap',[]);
for i = 1:size(acc,2)
    sens.acc = acc(:,i);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.4 0.3 0.5]);
    ft_topoplotER(cfg,sens);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.4 0.3 0.5]);
    c = colorbar; c.Label.String = dec_args.result_type;
    if ~isempty(dec_args.time)
        text(-0.5,-0.5, pad(num2str(round(dec_args.time(i),3)),10), 'FontWeight', 'normal');
    end;
    F(i) = getframe(gcf);
end;

vid_obj = VideoWriter(output_file);
vid_obj.FrameRate=2;
open(vid_obj);
writeVideo(vid_obj,F)
close
end

