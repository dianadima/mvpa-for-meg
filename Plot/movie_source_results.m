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
%   'style', default 'searchlight' (can be 'searchlight' or 'centroid')
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
addParameter(p, 'view', [0 90]); %view - default = from above
addParameter(p, 'result_type', 'Accuracy (%)');
addParameter(p, 'framerate', 2);
addParameter(p, 'style', 'searchlight');
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
sourcemodel.coordsys = 'mni';
if isempty(source_idx)
    
    if size(sourcemodel.pos(sourcemodel.inside,:,:),1) ~= size(acc,1)
        error('Please provide source indices or ensure accuracy dimension 1 fits number of inside sources in FT sourcemodel.')
    end;
    
    pow = NaN(size(sourcemodel.pos,1),size(acc,2));
    pow(sourcemodel.inside,:) = acc;
    
else
   
    all_acc = nan(length(find(sourcemodel.inside==1)),size(acc,2));
    
    for t = 1:size(acc,2)
        for i = 1:length(source_idx)
            if ~isnan(accuracy(i,t)) && accuracy(i,t)~=0
                if length(source_idx)==size(all_acc,1) && strcmp(p.Results.style, 'centroid')
                    all_acc(i,t) = acc(i,t);
                else
                    all_acc(source_idx{i},t) = acc(i,t);
                end
            end
        end
    end
    
    pow = NaN(size(sourcemodel.pos,1),size(acc,2));
    pow(sourcemodel.inside,:) = all_acc;
    
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


F(size(acc,2)) = struct('cdata',[],'colormap',[]);
for i = 1:size(accuracy,2)
    sourcemodel.pow = pow(:,i);
    if ~isempty(p.Results.roi), cfg.roi = p.Results.roi{i}; end;
    ft_sourceplot(cfg,sourcemodel); view(p.Results.view);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0.4 0.4 0.6]);
    c = colorbar; c.Label.String = p.Results.result_type;
    F(i) = getframe(gcf);
    close;    
end;

vid_obj = VideoWriter(output_file);
vid_obj.FrameRate = p.Results.framerate;
open(vid_obj);
writeVideo(vid_obj,F)

end

