function [] = get_sensor_info( dataset, output_file, varargin )
% Inputs: dataset, output file path. Optional: 
%  decoding_window: decoding window or window of
%   interest (to save time axis for use in plotting etc.). If empty, whole
%   axis saved.
%  sample_rate: sampling rate (consider whether you'll be resampling in your analysis;
%   this should be the final sample rate; default 1200).
% Output: sensor neighbour structure + channel ids + time axis, saved to
% .mat file that will be used in further decoding scripts.
% Uses fieldtrip toolbox (Oostenveld et al., 2011).
% DCD 2017

p = inputParser;
addParameter(p, 'decoding_window', []);
addParameter(p, 'sample_rate', 1200);
parse(p, varargin{:});

hdr = ft_read_header(dataset);

%get neighbour structure for configuration
cfg = [];
cfg.method = 'template';
cfg.template = [hdr.grad.type '_neighb.mat'];
cfg.layout = [hdr.grad.type '.lay'];
neighbours = ft_prepare_neighbours(cfg, hdr.grad); %#ok<NASGU> %neighblabel contains the neighbours of every single channel

%select MEG channels
for i = 1:length(hdr.label)
    idx = cellfun('isempty',strfind(hdr.label, 'M'));
end;
chan_labels = hdr.label(~idx);

if isempty(p.Results.decoding_window)
    time = -(hdr.nSamplesPre/p.Results.sample_rate):(1/p.Results.sample_rate):((hdr.nSamples-hdr.nSamplesPre)/p.Results.sample_rate);
else
    time = p.Results.decoding_window:(1/p.Results.sample_rate):p.Results.decoding_window(2);
end;

fprintf('\nSaving neighbourhood structure and labels for %d channels...', length(chan_labels));
fprintf('\nSaving time axis from %f to %f...\n', time(1), time(end));
save(output_file, 'neighbours', 'chan_labels', 'time'); %all useful things

end

