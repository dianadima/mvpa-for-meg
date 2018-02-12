function [ ] = plot_source_space( accuracy, source_idx, varargin )
% Plot source-space decoding results.
% Inputs: accuracy (or other metric of choice); must be a vector with the right number of sources.
%         source_idx: cell array obtained from get_source_info, shows where
%         each source belongs in a searchlight approach. Can be left empty
%         if your accuracy vector contains all sources.
% Optional inputs:
%         'sourcemodel', default is the template 10 mm Fieldtrip sourcemodel. Can specify any full path to a sourcemodel file.
%         'colormap', default 'jet'.
%         'colorlim', default [40 100].
%         'visible', default 'on' - open figure window or not. 

[~, ftdir] = ft_version; %get FT directory

p = inputParser;
addParameter(p, 'sourcemodel', fullfile(ftdir, 'template', 'sourcemodel', 'standard_sourcemodel10mm'));
addParameter(p, 'colormap', 'jet');
addParameter(p, 'colorlim', [40 100]);
addParameter(p, 'visible', 'on');
parse(p, varargin{:});

if size(accuracy,1)>1
    accuracy = accuracy';
end;

%load sourcemodel
if ischar(p.Results.sourcemodel)
    [~,~,ext] = fileparts(p.Results.sourcemodel);
    if strcmp(ext, '.mat')
        load(p.Results.sourcemodel);
    else
        sourcemodel = ft_read_headshape(p.Results.sourcemodel);
    end;
else
    sourcemodel = p.Results.sourcemodel;
end;
if ~isfield(sourcemodel, 'inside')
    sourcemodel.inside = true(size(sourcemodel.pos,1),1);
end;
sourcemodel = ft_convert_units(sourcemodel, 'mm');

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
cfg.visible = p.Results.visible;

ft_sourceplot(cfg,sourcemodel)


end

