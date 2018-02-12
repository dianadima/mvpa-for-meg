classdef decoding_args
    
    properties (Access = public)
        
        %default settings for time-resolved decoding
        channels = 'MEG';
        decoding_window = [];
        window_length = 1;
        time = [];
        
    end
    
    methods
        
        function obj = decoding_args(varargin)
            
            if ~isempty(varargin)
                obj.channels = varargin{1};
                obj.decoding_window = varargin{2};
                obj.window_length = varargin{3};
                obj.time = varargin{4};
            end;
        end
        
    end
        
end
    