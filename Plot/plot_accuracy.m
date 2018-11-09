function [ ] = plot_accuracy( accuracy, errorbar, varargin )
% Plot decoding accuracy (or other metric) over time, using transparent shading for the error bar.
% Function can be called in a loop to plot multiple accuracies on top of each other.
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

%get or set time axis
if isempty(p.Results.time)
    time = 1:length(accuracy);
else
    time = p.Results.time;
end;

if length(time)~=size(errorbar,2) || length(accuracy)~=size(errorbar,2)
    error('Time, accuracy and error should have the same length');
end;

%get error bars
if size(errorbar,1)==1
        upper = accuracy + errorbar;
        lower = accuracy - errorbar;
elseif size(errorbar,1)==2
    if errorbar(1,1)>errorbar(2,1)
         errorbar = flipud(errorbar);
    end;
    upper = errorbar(2,:);
    lower = errorbar(1,:);
else
    error('The errorbar can contain one or two row vectors');
end;

%plot time-course with shaded error bar
l0 = patch([time fliplr(time)], [smooth(upper, p.Results.smooth)' fliplr(smooth(lower, p.Results.smooth)')], p.Results.color, 'EdgeColor', 'none'); hasbehavior(l0, 'legend', false);
alpha 0.15; hold on; 
if ~isempty(p.Results.chance)
    l1 = line([time(1) time(end)], [p.Results.chance p.Results.chance], 'color', [0.5 0.5 0.5]); hold on; hasbehavior(l1, 'legend', false);
end
if ~isempty(time==0)
    l2 = line([0 0], [p.Results.ylim(1) p.Results.ylim(2)], 'color', [0.5 0.5 0.5]); hasbehavior(l2, 'legend', false);
end;
if ~isempty(p.Results.legend)
    pl = plot(time, smooth(accuracy, p.Results.smooth), 'color',p.Results.color, 'LineWidth', p.Results.linewidth,'LineStyle', p.Results.linestyle, 'DisplayName', p.Results.legend); hold on; hasbehavior(pl, 'legend', true);
else
    plot(time, smooth(accuracy, p.Results.smooth), 'color',p.Results.color, 'LineWidth', p.Results.linewidth,'LineStyle', p.Results.linestyle); hold on; 
end;

%mark significant time points with horizontal line & vertical line for onset
if ~isempty(p.Results.signif)
    l3 = plot(p.Results.signif, p.Results.signif_ylocation, 'Marker', 's', 'MarkerSize', p.Results.signif_marker_size, 'MarkerFaceColor', p.Results.color,'MarkerEdgeColor', p.Results.color);
    for j = 1:length(l3), hasbehavior(l3(j), 'legend', false); end; hold on;
    l4 = line([p.Results.signif(1,1) p.Results.signif(1,1)], p.Results.ylim, 'color', p.Results.color); hold on; hasbehavior(l4, 'legend', false); hold on;
end;

ylim(p.Results.ylim);
if isempty(p.Results.xlim)
    xlim([time(1) time(end)]);
else
    xlim([p.Results.xlim(1) p.Results.xlim(2)]);
end;
ylabel(p.Results.ylabel);
xlabel(p.Results.xlabel);
box off;

if ~isempty(p.Results.legend)
    legend('-DynamicLegend');
    legend boxoff;
end;

set(gca,'FontSize',12);

end

