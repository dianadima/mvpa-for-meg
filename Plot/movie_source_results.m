function [ vid_obj ] = movie_source_results( results, output_file, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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
    error('Results should be a 2d or 3d matrix containing subjects x time x channels');
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

F(size(acc,1)) = struct('cdata',[],'colormap',[]);
for i = 1:size(results,1)
    sourcemodel.acc(sourcemodel.inside) = acc(i,:)';
    ft_sourceplot(cfg,sourcemodel);
    F(i) = getframe(gca);
end;

vid_obj = VideoWriter(output_file);
vid_obj.FrameRate=2;
open(vid_obj);
writeVideo(vid_obj,F)

end

