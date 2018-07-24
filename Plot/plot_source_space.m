function [ ] = plot_source_space( accuracy, varargin )
% Plot source-space decoding results.
% Inputs: accuracy (or other metric of choice); must be a vector with the right number of sources.
%         
% Optional inputs:
%         'source_idx': cell array obtained from get_source_info, shows where
%                       each source belongs in a searchlight/AAL approach.
%         'sourcemodel', default is the template 10 mm Fieldtrip sourcemodel. Can specify any full path to a sourcemodel file.
%         'colormap', default 'jet'.
%         'colorlim', default [40 100].
%         'visible', default 'on' - open figure window or not. 
%         'roi', default [], mask all ROIs except this (AAL ROI string)
%         'inflated', default true - inflated surface or regular
%         'hemisphere', default 'both' (can be 'right', 'left', or 'both')
%         'style', default 'searchlight' (can be 'searchlight' or 'centroid')

[~, ftdir] = ft_version; %get FT directory

p = inputParser;
addParameter(p, 'source_idx', []);
addParameter(p, 'sourcemodel', fullfile(ftdir, 'template', 'sourcemodel', 'standard_sourcemodel3d10mm.mat'));
addParameter(p, 'colormap', 'jet');
addParameter(p, 'colorlim', [40 100]);
addParameter(p, 'visible', 'on');
addParameter(p, 'surface_file', []);
addParameter(p, 'inflated', true);
addParameter(p, 'roi', []);
addParameter(p, 'hemisphere', 'both');
addParameter(p, 'view', [0 90]);
addParameter(p, 'style', 'searchlight');
parse(p, varargin{:});

if size(accuracy,1)>1 && size(accuracy,2)>1
    error('The parameter you wish to plot must be a vector')
end
    
%load sourcemodel
if ischar(p.Results.sourcemodel)
    [~,~,ext] = fileparts(p.Results.sourcemodel);
    if strcmp(ext, '.mat')
        load(p.Results.sourcemodel);
    else
        sourcemodel = ft_read_headshape(p.Results.sourcemodel);
    end
else
    sourcemodel = p.Results.sourcemodel;
end
if ~isfield(sourcemodel, 'inside')
    sourcemodel.inside = true(size(sourcemodel.pos,1),1);
end
sourcemodel = ft_convert_units(sourcemodel, 'mm');
sourcemodel.coordsys = 'mni';

source_idx = p.Results.source_idx;

if isempty(source_idx)
    
    if size(sourcemodel.pos(sourcemodel.inside,:,:),1) ~= length(accuracy)
        error('Please provide source indices or ensure accuracy vector length fits number of inside sources in FT sourcemodel.')
    end
    
    sourcemodel.pow = NaN(1,size(sourcemodel.pos,1));
    sourcemodel.pow(sourcemodel.inside) = accuracy;
    
else
    
    pow = nan(1,length(find(sourcemodel.inside==1)));
    for i = 1:length(pow)
        if ~isnan(accuracy(i)) && accuracy(i)~=0
            if strcmp(p.Results.style, 'centroid')
                    pow(i) = accuracy(i); %assign to centroids
            elseif strcmp(p.Results.style, 'searchlight')
                     pow(source_idx{i}) = accuracy(i);
            end
        end
    end
   
    sourcemodel.pow = nan(length(sourcemodel.inside),1);
    sourcemodel.pow(sourcemodel.inside==1) = pow;
    sourcemodel.coordsys = 'mni';
    
end
    
cfg = [];
cfg.method = 'surface'; 
cfg.funparameter = 'pow';
if isempty(p.Results.roi)
    cfg.maskparameter = 'pow';
end
cfg.funcolormap = p.Results.colormap;
cfg.funcolorlim = p.Results.colorlim;
cfg.opacitylim = p.Results.colorlim;
cfg.opacitymap = 'rampup';
cfg.projmethod = 'nearest';
cfg.projcomb = 'mean';
cfg.projvec = 1;
cfg.projweight = 1;
if ~isempty(p.Results.roi)
    cfg.atlas = ft_read_atlas([ftdir '/template/atlas/aal/ROI_MNI_V4.nii']);
    cfg.roi = p.Results.roi;
end
if isempty(p.Results.surface_file)
    cfg.surffile = ['surface_white_' p.Results.hemisphere '.mat'];
    if p.Results.inflated
        cfg.surfinflated = ['surface_inflated_' p.Results.hemisphere '.mat'];
    end
else
    cfg.surffile = p.Results.surface_file;
end
cfg.camlight = 'no';
cfg.visible = p.Results.visible;

ft_sourceplot(cfg,sourcemodel); view(p.Results.view);
if ~p.Results.inflated,  camlight('headlight','infinite'); material('dull'); end;

end

