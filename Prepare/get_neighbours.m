function [] = get_neighbours( data, configuration, output_file )
% Inputs: data, MEG sensor configuration, output file path.
% Output: sensor neighbour structure + channel ids + time axis, saved to
% .mat file that will be used in further decoding scripts.
% Uses fieldtrip toolbox (Oostenveld et al., 2011).
% DCD 2017

cfg = [];
cfg.method = 'template';
cfg.template = [configuration '_neighb.mat'];
cfg.layout = [configuration '.lay'];
neighbours = ft_prepare_neighbours(cfg, data); %neighblabel contains the neighbours of every single channel

chan_labels = data.label; time = data.time{1};
save(output_file,'neighbours', 'chan_labels', 'time'); %all useful things

end

