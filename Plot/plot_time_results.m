function [ ] = plot_time_results( accuracy, errorbar, varargin )
% Plot decoding accuracy (or other metric) over time, using transparent shading for the error bapr.
% Function can be called in a loop to plot multiple accuracies on top of each othepr.
%
% Inputs: accuracy and errorbar (e.g., standard deviation or SEM to be added to & subtracted from accuracy). Both same length
% Optional inputs:
%       'time', time axis; will be plotted on X axis and needs to match length of accuracy and errorbar vectors.
%       'signif', significant time points, will be underlined with horizontal lines and onset marked with vertical bars. Must be in same unit as x-axis
%       'signif_ylocation', location of vertical lines marking significance on the y-axis (default 47). Especially useful if calling this function in a loop.
%       'signif_marker_size', controls the size of the line/markers denoting significance (default 4).
%       'color', accuracy & errorbar color, default 'k', black.
%       'linewidth', accuracy line width (default 1.5)
%       'smooth', factor by which to smooth the accuracy with a moving average (default 1 - no smoothing)
%       'ylim', default [40 100]
%       'xlim', default [] (data-driven)
%       'ylabel', default 'Accuracy (%)'
%       'xlabel', default 'Time (s)'
%       'legend' - specify string corresponding to the accuracy plotted. If function is called in a loop, 
%                  legends will be appended, such that the final plot will have the complete legend of plotted accuracies. Default, [] (none). (Semidocumented function DynamicLegend used here)    
%
% DC Dima 2018 (diana.c.dima@gmail.com)

%parse inputs and set defaults for name-value pairs
p = inputParser;
addParameter(p, 'time',[]);
addParameter(p, 'color', 'k');
addParameter(p, 'signif', []);
addParameter(p, 'signif_ylocation', 47);
addParameter(p, 'signif_marker_size',4);
addParameter(p, 'smooth', 1);
addParameter(p, 'ylim', [40 100]);
addParameter(p, 'chance', 50);
addParameter(p, 'xlim', []);
addParameter(p, 'linewidth', 1.5);
addParameter(p, 'linestyle', '-');
addParameter(p, 'ylabel', 'Accuracy (%)');
addParameter(p, 'xlabel', 'Time (s)');
addParameter(p, 'legend', []);
parse(p, varargin{:});

pr = p.Results;

%get or set time axis
if isempty(pr.time)
    time = 1:length(accuracy);
else
    time = pr.time;
end

%column vectors
if size(time,1)~=1, time = time'; end
if size(errorbar,1)~=1, errorbar = errorbar'; end
if size(accuracy,1)~=1, accuracy = accuracy'; end

if size(time,2)~=size(errorbar,2) || size(accuracy,2)~=size(errorbar,2)
    error('Time, accuracy and error should have the same length');
end

%get error bars
if size(errorbar,1)==1
        upper = accuracy + errorbar;
        lower = accuracy - errorbar;
elseif size(errorbar,1)==2
    if errorbar(1,1)>errorbar(2,1)
         errorbar = flipud(errorbar);
    end
    upper = errorbar(2,:);
    lower = errorbar(1,:);
else
    error('The errorbar can contain one or two row vectors');
end

%plot time-course with shaded error bar
l0 = patch([time fliplr(time)], [smooth(upper, pr.smooth)' fliplr(smooth(lower, pr.smooth)')], pr.color, 'EdgeColor', 'none'); hasbehavior(l0, 'legend', false);
alpha 0.15; hold on; 
if isempty(pr.ylim), pr.ylim = get(gca, 'ylim'); end
if ~isempty(pr.chance)
    l1 = line([time(1) time(end)], [pr.chance pr.chance], 'color', [0.5 0.5 0.5]); hold on; hasbehavior(l1, 'legend', false);
end
if ~isempty(time==0)
    l2 = line([0 0], [pr.ylim(1) pr.ylim(2)], 'color', [0.5 0.5 0.5]); hasbehavior(l2, 'legend', false);
end
if ~isempty(pr.legend)
    pl = plot(time, smooth(accuracy, pr.smooth), 'color',pr.color, 'LineWidth', pr.linewidth,'LineStyle', pr.linestyle, 'DisplayName', pr.legend); hold on; hasbehavior(pl, 'legend', true);
else
    plot(time, smooth(accuracy, pr.smooth), 'color',pr.color, 'LineWidth', pr.linewidth,'LineStyle', pr.linestyle); hold on; 
end

%mark significant time points with horizontal line & vertical line for onset
if ~isempty(pr.signif) && sum(pr.signif)~=0
    if length(pr.signif)==length(time) && max(pr.signif)==1, pr.signif = time(logical(pr.signif)); end %logical indexing case
    l3 = plot(pr.signif, pr.signif_ylocation, 'Marker', 's', 'MarkerSize', pr.signif_marker_size, 'MarkerFaceColor', pr.color,'MarkerEdgeColor', pr.color);
    for j = 1:length(l3), hasbehavior(l3(j), 'legend', false); end; hold on;
    l4 = line([pr.signif(1,1) pr.signif(1,1)], pr.ylim, 'color', pr.color); hold on; hasbehavior(l4, 'legend', false); hold on;
end

if ~isempty(pr.ylim), ylim(pr.ylim); end
if isempty(pr.xlim)
    xlim([time(1) time(end)]);
else
    xlim([pr.xlim(1) pr.xlim(2)]);
end
ylabel(pr.ylabel);
xlabel(pr.xlabel);
box off

if ~isempty(pr.legend)
    legend('-DynamicLegend')
    legend boxoff
end

set(gca,'FontSize',14);

end

