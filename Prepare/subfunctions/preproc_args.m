classdef preproc_args
    
    properties (Access = public)
        
        %define default values for function parameters
        prestim = 0.99;
        poststim = 0.99;
        baseline = [-0.5 0];
        bandpass_filter = [0 100];
        resamplefs = 600;
        
    end
    
    methods
        
        function obj = preproc_args(varargin)
            
            if ~isempty(varargin)
                obj.prestim = varargin{1};
                obj.poststim = varargin{2};
                obj.baseline = varargin{3};
                obj.bandpass_filter = varargin{4};
                obj.resamplefs = varargin{5};
            end;
                
        end;
    end;
    
end
        

