function plot_temporal_generalization(results, varargin)
% Plots results of temporal generalization decoding.
% Inputs: results: matrix of accuracy/decoding performance. Must be time x time, or subjects x time x time.
%
% Optional inputs:
%   'time'(default []), time axis, if you wish to specify time units, rather than sampled points.
%   'mask', binary mask (significant = 1). Significant clusters will be marked with white contours  
%   'clustersize' (default 1), minimal cluster size in the binary mask (2D) to be plotted
%   'colorlim' (default [40 100]): colour limits
%   'colormap' (default 'jet')
%   'colorbar_label' (default 'Accuracy (%)')
%   'title' (default []), figure title
%
% DC Dima 2019 (diana.c.dima@gmail.com)

p = inputParser;
addParameter(p, 'colorlim', [40 100]);
addParameter(p, 'colormap', 'jet');
addParameter(p, 'time', []);
addParameter(p, 'mask', []);
addParameter(p, 'clustersize',1)
addParameter(p, 'clusters', []);
addParameter(p, 'title', []);
addParameter(p, 'colorbar', true);
addParameter(p, 'colorbar_label', 'Accuracy (%)');
addParameter(p, 'grid', []);
parse(p, varargin{:});
opt = p.Results;

%create time axis
if ~isempty(opt.time)
    time = opt.time;
else
    if ismatrix(results)
        time = 1:size(results,2);
    elseif ndims(results)==3
        time = 1:size(results,3);
    end;
end;

if ismatrix(results)
    acc = results;
elseif ndims(results)==3
    acc = squeeze(mean(results,1));
    fprintf('Warning: assuming subjects are 1st dimension of accuracy matrix....')
else
    error('Results should be a 2d or 3d matrix containing subjects x time x time');
end;

%have 6 ticklabels for readability
nticks = 6;
ntp = length(time);
nstep = floor(ntp/nticks);

imagesc(acc)
colormap(opt.colormap)
if ~isempty(opt.colorlim), caxis(opt.colorlim); end
f = gca;
set(f, 'xtick', 1:nstep:ntp)
set(f,'xticklabels', round(time(1:nstep:ntp),1))
set(f,'ytick', 1:nstep:ntp)
set(f,'yticklabels', round(time(1:nstep:ntp),1))
set(f, 'ydir', 'normal')
if ~isempty(opt.title), title(opt.title, 'FontWeight', 'normal'); end
set(f, 'FontSize', 14)
xlabel('Test time')
ylabel('Train time')
if opt.colorbar
    c = colorbar; c.Location = 'eastoutside';
    c.Label.String = opt.colorbar_label;
end
hold on

clusters = [];
if ~isempty(opt.mask)
    comp = bwconncomp(opt.mask); %get the clusters
    clusters = comp.PixelIdxList;
elseif ~isempty(opt.clusters)
    clusters = opt.clusters;
end

if ~isempty(clusters)
    for j = 1:length(clusters)
        if length(clusters{j})>=opt.clustersize
            z = zeros(size(acc)); 
            z(clusters{j}) = 1;
            [yy,xx] = find(z==1);
            k = boundary(xx,yy,1);
            plot(xx(k), yy(k), 'w', 'LineWidth', 1.5);
            hold on
        end
    end
end

if ~isempty(opt.grid)
    grid on;
    ax = gca;
    ax.GridColor = opt.grid;
end


end

