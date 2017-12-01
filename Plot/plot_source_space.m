function [ ] = plot_source_space( accuracy, source_idx, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
addParameter(p, 'sourcemodel_resolution',10);
addParameter(p, 'colormap', 'jet');
addParameter(p, 'colorlim', [40 100]);
parse(p, varargin{:});

if size(accuracy,1)>1
    accuracy = accuracy';
end;

%load sourcemodel
[~, ftdir] = ft_version; %get FT directory
if p.Results.sourcemodel_resolution==7.5
    load(fullfile(ftdir, 'template', 'sourcemodel', 'standard_sourcemodel3d7point5mm'));
else
    load(fullfile(ftdir, 'template', 'sourcemodel', ['standard_sourcemodel3d' num2str(p.Results.sourcemodel_resolution) 'mm']));
end;
sourcemodel = ft_convert_units(sourcemodel, 'mm'); %#ok<NODEF>

if isempty(source_idx)
    
    if size(sourcemodel.pos(sourcemodel.inside,:,:),1) ~= length(accuracy)
        error('Please provide source indices or ensure accuracy vector length fits number of inside sources in FT sourcemodel.')
    end;
    
    sourcemodel.pow = NaN(1,size(sourcemodel.pos,1));
    sourcemodel.pow(sourcemodel.inside) = accuracy;
    
else
    
    idx = unique(cell2mat(source_idx));
    inside = 1:size(sourcemodel.pos(sourcemodel.inside,:,:),1);
    outside_idx = inside; outside_idx(idx) = [];
    
    for i = 1:length(source_idx)
        inside(ismember(inside,source_idx{i})) = accuracy(i);
    end;
    
    inside(outside_idx) = NaN;
    
    sourcemodel.pow = NaN(1,size(sourcemodel.pos,1));
    sourcemodel.pow(sourcemodel.inside) = inside;
    
end;
    
cfg = [];
cfg.method = 'surface'; 
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap = p.Results.colormap;
cfg.funcolorlim = p.Results.colorlim;
cfg.opacitylim = p.Results.colorlim;
cfg.opacitymap = 'rampup';
cfg.projmethod = 'nearest';
cfg.surffile = 'surface_white_both.mat';
cfg.camlight = 'no';

ft_sourceplot(cfg,sourcemodel)


end

