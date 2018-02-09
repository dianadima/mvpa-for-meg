function [ ] = plot_accuracy( accuracy, errorbar, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p = inputParser;
addParameter(p, 'time',[]);
addParameter(p, 'color', 'k');
addParameter(p, 'signif', []);
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

l0 = patch([time fliplr(p.Results.time)], [smooth(upper, p.Results.smooth)' fliplr(smooth(lower, p.Results.smooth)')], p.Results.color, 'EdgeColor', 'none'); hasbehavior(l0, 'legend', false);
alpha 0.15; hold on; 
l1 = line([time(1) p.Results.time(end)], [50 50], 'color', [0.5 0.5 0.5]); hold on; hasbehavior(l1, 'legend', false);
if ~isempty(time==0)
    l2 = line([0 0], [p.Results.ylim(1) p.Results.ylim(2)], 'color', [0.5 0.5 0.5]); hasbehavior(l2, 'legend', false);
end;
pl = plot(time, smooth(accuracy, p.Results.smooth), 'color',p.Results.color, 'LineWidth', p.Results.linewidth, 'DisplayName', p.Results.legend); hold on; hasbehavior(pl, 'legend', true);

if ~isempty(p.Results.signif)
    for j = 1:size(p.Results.signif,1)
        l3 = line([p.Results.signif(j,1) p.Results.signif(j,2)], p.Results.ylim+5, 'color', p.Results.color, 'LineWidth', 2); hold on; hasbehavior(l3, 'legend', false);
    end;
    l4 = line([p.Results.signif(1,1) p.Results.signif(1,1)], p.Results.ylim, 'color', p.Results.color); hold on; hasbehavior(l4, 'legend', false);
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

