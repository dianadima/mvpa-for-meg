classdef decoding_args
    
    properties (Access = public)
        
        %default settings for time-resolved decoding
        channels = 'MEG';
        decoding_window = [-0.2 0.9];
        window_length = 1;
        
    end
    
    methods
        
        function obj = decoding_args(varargin)
            
            if ~isempty(varargin)
                obj.channels = varargin{1};
                obj.decoding_window = varargin{2};
                obj.window_length = varargin{3};
            end;
        end
        
    end
        
end
    