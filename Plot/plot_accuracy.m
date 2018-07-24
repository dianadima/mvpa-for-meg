function [ ] = plot_accuracy( accuracy, errorbar, varargin )
% Plot decoding accuracy (or other metric) over time, using transparent shading for the error bar.
% Function can be called in a loop to plot multiple accuracies on top of each other.
%
% Inputs: accuracy and errorbar (e.g., standard deviation or SEM to be added to & subtracted from accuracy). Both same length
% Optional inputs:
%       'time', time axis; will be plotted on X axis and needs to match length of accuracy and errorbar vectors.
%       'signif', significant time points, will be underlined with horizontal lines and onset marked with vertical bars. Must be in same unit as x-axis
%       'signif_ylocation', location of vertical lines marking significance on the y-axis (default 47). Especially useful if calling this function in a loop.
%       'color', accuracy & errorbar color, default 'k', black.
%       'linewidth', accuracy line width (default 1.5)
%       'smooth', factor by which to smooth the accuracy with a moving average (default 1 - no smoothing)
%       'ylim', default [40 100]
%       'xlim', default [] (data-driven)
%       'ylabel', default 'Accuracy (%)'
%       'xlabel', default 'Time (s)'
%       'legend' - specify string corresponding to the accuracy plotted. If function is called in a loop, 
%                   legends will be appended, such that the final plot will have the complete legend of plotted accuracies. Default, [] (none). (Semidocumented function DynamicLegend used here)    

p = inputParser;
addParameter(p, 'time',[]);
addParameter(p, 'color', 'k');
addParameter(p, 'signif', []);
addParameter(p,'signif_ylocation', 47);
addParameter(p, 'smooth', 1);
addParameter(p, 'ylim', [40 100]);
addParameter(p, 'xlim', []);
addParameter(p, 'linewidth', 1.5);
addParameter(p, 'ylabel', 'Accuracy (%)');
addParameter(p, 'xlabel', 'Time (s)');
addParameter(p, 'legend', []);
parse(p, varargin{:});

if isempty(p.Results.time)
    time = 1:length(accuracy);
else
    time = p.Results.time;
end;
    
if length(time)~=size(errorbar,2) || length(accuracy)~=size(errorbar,2)
    error('Time, accuracy and error should have the same length');
end;

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

l0 = patch([time fliplr(time)], [smooth(upper, p.Results.smooth)' fliplr(smooth(lower, p.Results.smooth)')], p.Results.color, 'EdgeColor', 'none'); hasbehavior(l0, 'legend', false);
alpha 0.15; hold on; 
l1 = line([time(1) time(end)], [50 50], 'color', [0.5 0.5 0.5]); hold on; hasbehavior(l1, 'legend', false);
if ~isempty(time==0)
    l2 = line([0 0], [p.Results.ylim(1) p.Results.ylim(2)], 'color', [0.5 0.5 0.5]); hasbehavior(l2, 'legend', false);
end;
if ~isempty(p.Results.legend)
    pl = plot(time, smooth(accuracy, p.Results.smooth), 'color',p.Results.color, 'LineWidth', p.Results.linewidth, 'DisplayName', p.Results.legend); hold on; hasbehavior(pl, 'legend', true);
else
    plot(time, smooth(accuracy, p.Results.smooth), 'color',p.Results.color, 'LineWidth', p.Results.linewidth); hold on; 
end;

if ~isempty(p.Results.signif)
    for j = 1:length(p.Results.signif)
        tmp = find(time==p.Results.signif(j));  
        if time(end)>p.Results.signif(j), tmp = time(tmp+1); end; 
        l3 = line([p.Results.signif(j) tmp], [p.Results.signif_ylocation p.Results.signif_ylocation], 'color', p.Results.color, 'LineWidth', 2); hasbehavior(l3, 'legend', false); hold on; 
    end;
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

