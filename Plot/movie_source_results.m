function [] = movie_source_results( results, output_file, varargin )
% Plot sensor-space searchlight decoding results as a movie.
% Inputs: results: matrix of accuracy/decoding performance. Must be channels x time, or subjects x channels x time.
%         neighbours: sensor grouping structure obtained using get_sensor_info (i.e., fieldtrip function prepare_neighbours).
%         output_file: movie filename
% Optional inputs:
%   'configuration' (default 'CTF275'): specify configuration of sensors. The fieldtrip layout <configuration>.lay will be loaded for plotting.
%   'clim' (default [40 100]): colour limits
%   'colormap' (default 'jet')
%   'result_type' (default 'Accuracy (%)'): will be plotted as colorbar axis
% Need to add capability for ROI and searchlight plotting

p = inputParser;
addParameter(p, 'sourcemodel_resolution',10);
parse(p, varargin{:});

[~, ftdir] = ft_version; %get FT directory

if p.Results.sourcemodel_resolution==7.5
    load(fullfile(ftdir, 'template', 'sourcemodel', 'standard_sourcemodel3d7point5mm'));
else
    load(fullfile(ftdir, 'template', 'sourcemodel', ['standard_sourcemodel3d' num2str(p.Results.sourcemodel_resolution) 'mm']));
end;

sourcemodel = ft_convert_units(sourcemodel, 'mm'); %#ok<NODEF>
sourcemodel.acc = NaN(1,size(sourcemodel.pos,1));

if ismatrix(results)
    acc = results;
elseif ndims(results)==3
    acc = squeeze(mean(results,1));
    fprintf('Warning: assuming subjects are 1st dimension of accuracy matrix....')
else
    error('Results should be a 2d or 3d matrix containing subjects x channels x time');
end;

cfg = [];
cfg.method = 'surface';
cfg.funparameter = 'acc';
cfg.maskparameter = 'acc';
cfg.funcolormap = 'jet';
cfg.opacitylim = 'maxabs';
cfg.opacitymap = 'vdown';
cfg.projmethod = 'nearest';
cfg.surffile = 'surface_white_both.mat';
cfg.visible = 'off';
cfg.camlight = 'no';

F(size(acc,2)) = struct('cdata',[],'colormap',[]);
for i = 1:size(results,2)
    sourcemodel.acc(sourcemodel.inside) = acc(:,i);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.4 0.3 0.5]);
    ft_sourceplot(cfg,sourcemodel);
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

end

