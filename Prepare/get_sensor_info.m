function [neighbours] = get_sensor_info( dataset, varargin )
% Inputs: dataset, output file path. Optional: 
%  decoding_window: decoding window or window of
%         interest (to save time axis for use in plotting etc.). If empty, whole
%         axis saved.
%  sample_rate: sampling rate - If you are going to resample the data, this
%         need to be provided, otherwise you'll get a time axis based on the dataset sampling rate.
% Output: sensor neighbour structure (includes labels, neighbours, and time).
% Uses Fieldtrip toolbox (Oostenveld et al., 2011).

hdr = ft_read_header(dataset);

%get neighbour structure for configuration
cfg = [];
cfg.method = 'template';
cfg.template = [hdr.grad.type '_neighb.mat'];
cfg.layout = [hdr.grad.type '.lay'];
neighbours = ft_prepare_neighbours(cfg, hdr.grad);%neighblabel contains the neighbours of every single channel

p = inputParser;
addParameter(p, 'decoding_window', []);
addParameter(p, 'sample_rate', hdr.Fs);
parse(p, varargin{:});

if isempty(p.Results.decoding_window)
    [neighbours.time] = deal(-(hdr.nSamplesPre/hdr.Fs):(1/p.Results.sample_rate):((hdr.nSamples-hdr.nSamplesPre)/hdr.Fs));
else
    [neighbours.time] = deal(p.Results.decoding_window:(1/p.Results.sample_rate):p.Results.decoding_window(2));
end;


end

