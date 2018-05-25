function [] = rotating_source_results( accuracy, output_file, source_idx, varargin )
% Plot sensor-space searchlight decoding results as a movie.
% Inputs: results: matrix of accuracy/decoding performance. Must be channels x time, or subjects x channels x time.
%         neighbours: sensor grouping structure obtained using get_sensor_info (i.e., fieldtrip function prepare_neighbours).
%         output_file: movie filename
% Optional inputs:
%   'time' (default []) - time axis, to display next to each frame
%   'clim' (default [40 100]): colour limits
%   'colormap' (default 'jet')
%   'result_type' (default 'Accuracy (%)'): will be plotted as colorbar axis
% Need to add capability for ROI and searchlight plotting
[~, ftdir] = ft_version; %get FT directory

p = inputParser;
addParameter(p, 'sourcemodel', fullfile(ftdir, 'template', 'sourcemodel', 'standard_sourcemodel3d10mm.mat'));
addParameter(p, 'colormap', 'jet');
addParameter(p, 'colorlim', [40 100]);
addParameter(p, 'visible', 'on');
addParameter(p, 'inflated', true);
addParameter(p, 'hemisphere', 'both');
addParameter(p, 'roi', []);
addParameter(p, 'result_type', 'Accuracy (%)');
addParameter(p, 'framerate', 10);
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
    
    idx = unique(cell2mat(source_idx));
    inside = 1:size(sourcemodel.pos(sourcemodel.inside,:,:),1);
    outside_idx = inside; outside_idx(idx) = [];
    
    for i = 1:length(source_idx)
        inside(ismember(inside,source_idx{i})) = accuracy(i);
    end;
    
    inside(outside_idx) = NaN;
    
    sourcemodel.pow = NaN(1,size(sourcemodel.pos,1))';
    sourcemodel.pow(sourcemodel.inside) = inside;
    
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
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.4 0.3 0.5]);
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

