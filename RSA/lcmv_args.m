classdef lcmv_args
    
    properties
        
        normalize = true;
        marker = 'onset';
        bandpass = [0.1 100];
        prestim = 0.5;
        poststim = 0.9;
        baseline = [-0.5 0];
        resamplefs = 600;
        plot = 0;
        
    end
    
    methods
        
        function obj = lcmv_args(varargin)
            
            if ~isempty(varargin)
                obj.normalize = varargin{1};
                obj.marker = varargin{2};
                obj.bandpass = varargin{3};
                obj.prestim = varargin{4};
                obj.poststim = varargin{5};
                obj.baseline = varargin{6};
                obj.resamplefs = varargin{7};
                obj.plot = varargin{8};
            end;
            
        end
        
    end
    
end
        