function [] = rotating_source_results( accuracy, source_idx, output_file, varargin )
% Plots source-space searchlight or ROI decoding results on a rotating brain.
% Inputs: accuracy: must be a vector with length = number of sources/ROIs/searchlights.
%         source_idx: source grouping structure obtained using get_source_info (has ROI or searchlight voxel indices).
%         output_file: movie filename, will be saved in current directory
%
% Optional inputs:
%   'sourcemodel' (default 10 mm grid) - sourcemodel used in source reconstruction & plotting
%   'colorlim' (default [40 100]): colour limits
%   'colormap' (default 'jet')
%   'result_type' (default 'Accuracy (%)'): will be plotted as colorbar axis
%   'style' (default 'centroid'), can be 'searchlight' or 'centroid' - value assigned to cluster of neighbouring sources or only centroid
%   'roi' (default []), region of interest (AAL label)
%   'framerate' (default 10)
%   'view' (default [0 90]), default from above
%   'hemisphere' (default 'both'), can be 'right', 'left'
%   'inflated' (default true), inflated brain surface
%
% DC Dima 2018 (diana.c.dima@gmail.com)

[~, ftdir] = ft_version; %get FT directory

p = inputParser;
addParameter(p, 'sourcemodel', fullfile(ftdir, 'template', 'sourcemodel', 'standard_sourcemodel3d10mm.mat'));
addParameter(p, 'colormap', 'jet');
addParameter(p, 'colorlim', [40 100]);
addParameter(p, 'inflated', true);
addParameter(p, 'hemisphere', 'both');
addParameter(p, 'roi', []);
addParameter(p, 'result_type', 'Accuracy (%)');
addParameter(p, 'framerate', 10);
addParameter(p, 'style', 'centroid');
parse(p, varargin{:});

if size(accuracy,1)>1 && size(accuracy,2)>1
    error('The parameter you wish to plot must be a vector')
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
sourcemodel.coordsys = 'mni';
if isempty(source_idx)
    
    if size(sourcemodel.pos(sourcemodel.inside,:,:),1) ~= size(accuracy,1)
        error('Please provide source indices or ensure accuracy dimension 1 fits number of inside sources in FT sourcemodel.')
    end;
    
    sourcemodel.pow = NaN(1,size(sourcemodel.pos,1));
    sourcemodel.pow(sourcemodel.inside) = accuracy;
    
else
    
    pow = nan(1,length(find(sourcemodel.inside==1)));
    for i = 1:length(source_idx)
        if ~isnan(accuracy(i)) && accuracy(i)~=0
            if length(source_idx)==length(pow) && strcmp(p.Results.style, 'centroid')
                pow(i) = accuracy(i); %assign to centroids
            else
                pow(source_idx{i}) = accuracy(i);
            end
        end
    end
   
    sourcemodel.pow = nan(length(sourcemodel.inside),1);
    sourcemodel.pow(sourcemodel.inside) = pow;
    
end;

cfg = [];
cfg.method = 'surface'; 
cfg.funparameter = 'pow';
if isempty(p.Results.roi)
    cfg.maskparameter = 'pow';
end;
cfg.funcolormap = p.Results.colormap;
cfg.funcolorlim = p.Results.colorlim;
cfg.opacitylim = p.Results.colorlim;
cfg.opacitymap = 'rampup';
cfg.projmethod = 'project';
cfg.projvec = 3;
if ~isempty(p.Results.roi)
    cfg.atlas = ft_read_atlas([ftdir '/template/atlas/aal/ROI_MNI_V4.nii']);
end;
cfg.surffile = ['surface_white_' p.Results.hemisphere '.mat'];
if p.Results.inflated
    cfg.surfinflated = ['surface_inflated_' p.Results.hemisphere '.mat'];
end;
cfg.camlight = 'no';
cfg.visible = 'on';
if ~isempty(p.Results.roi)
    cfg.roi = p.Results.roi;
end;

views = repmat(0:5:360,1,3);
F(length(views)) = struct('cdata',[],'colormap',[]);
ft_sourceplot(cfg,sourcemodel); 
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.4 0.4 0.6]);
c = colorbar; c.Label.String = p.Results.result_type; 
for i = 1:length(views)
    view([views(i) 0]);
    F(i) = getframe(gcf);
end;

vid_obj = VideoWriter(output_file);
vid_obj.FrameRate = p.Results.framerate;
open(vid_obj);
writeVideo(vid_obj,F)

end

