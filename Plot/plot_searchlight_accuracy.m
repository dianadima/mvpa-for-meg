function [ ] = plot_searchlight_accuracy(accuracy, neighbours, varargin)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    dec_args = decoding_args;
    dec_args.window_length = 0.05; %just this one different than default
elseif nargin>=3
    dec_args.channels = varargin{1};
    if nargin>=4
        dec_args.decoding_window = varargin{2};
        if nargin==5
            dec_args.window_length = varargin{3};
        end;
    end;
else
    error('Unexpected number of arguments. Check help file!')
end;

load(neighbours);

if ~isempty(dec_args.decoding_window)
    lims = [find(round(time,3)==dec_args.decoding_window(1)) find(round(time,3)==dec_args.decoding_window(2))];
else
    lims = [1 length(time)];
end;

if ismatrix(accuracy)
    acc = accuracy;
elseif ndims(accuracy)==3
    acc = squeeze(mean(accuracy,1));
    fprintf('Warning: assuming subjects are 1st dimension of accuracy matrix....')
else
    error('Accuracy should be a 2d or 3d matrix containing subjects x time x channels');
end;

%create structure for plotting
sens.acc = acc;
sens.label = chan_labels;
sens.dimord = 'chan_time';
sens.time = time(lims(1):lims(2));
cfg.xlim = time(lims(1)):time_win:time(lims(2));
%cfg.zlim = [0.4 1];
cfg.layout = 'CTF275.lay';
cfg.parameter = 'acc';
cfg.style = 'straight';
ft_topoplotER(cfg,sens)

end

