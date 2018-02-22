function [] = movie_source_results( accuracy, output_file, source_idx, varargin )
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
addParameter(p, 'result_type', 'Accuracy (%)');
parse(p, varargin{:});

if ismatrix(accuracy)
    acc = accuracy;
elseif ndims(accuracy)==3
    acc = squeeze(mean(accuracy,1));
    fprintf('Warning: assuming subjects are 1st dimension of accuracy matrix....')
else
    error('Results should be a 2d or 3d matrix containing subjects x channels x time');
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
    
    if size(sourcemodel.pos(sourcemodel.inside,:,:),1) ~= size(accuracy,1)
        error('Please provide source indices or ensure accuracy dimension 1 fits number of inside sources in FT sourcemodel.')
    end;
    
    pow = NaN(size(sourcemodel.pos,1),size(accuracy,2));
    pow(sourcemodel.inside,:) = accuracy;
    
else
    
    idx = unique(cell2mat(source_idx));
    inside = 1:size(sourcemodel.pos(sourcemodel.inside,:,:),1);
    all_acc = repmat(inside,size(accuracy,2),1)';
    outside_idx = inside; outside_idx(idx) = [];
    
    for t = 1:size(accuracy,2)
        for i = 1:length(source_idx)
            inside(ismember(inside,source_idx{i})) = accuracy(i,t);
        end;
        all_acc(:,t) = inside';
        inside = 1:size(sourcemodel.pos(sourcemodel.inside,:,:),1);
    end;
    
    all_acc(outside_idx,:) = NaN;
    
    pow = NaN(size(sourcemodel.pos,1),size(accuracy,2));
    pow(sourcemodel.inside,:) = all_acc;
    
end;

cfg = [];
cfg.method = 'surface'; 
cfg.funparameter = 'pow';
cfg.maskparameter = 'pow';
cfg.funcolormap = p.Results.colormap;
cfg.funcolorlim = p.Results.colorlim;
cfg.opacitylim = p.Results.colorlim;
cfg.opacitymap = 'rampup';
cfg.projmethod = 'project';
cfg.projvec = 3;
cfg.surffile = ['surface_white_' p.Results.hemisphere '.mat'];
if p.Results.inflated
    cfg.surfinflated = ['surface_inflated_' p.Results.hemisphere '.mat'];
end;
cfg.camlight = 'no';
cfg.visible = 'on';


F(size(acc,2)) = struct('cdata',[],'colormap',[]);
for i = 1:size(accuracy,2)
    sourcemodel.pow = pow(:,i);
    ft_sourceplot(cfg,sourcemodel);
    c = colorbar; c.Label.String = p.Results.result_type;
    F(i) = getframe(gcf);
end;

vid_obj = VideoWriter(output_file);
vid_obj.FrameRate = 2;
open(vid_obj);
writeVideo(vid_obj,F)

end

