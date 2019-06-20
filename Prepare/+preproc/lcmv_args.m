% class containing LCMV beamforming & preprocessing parameters
% DC Dima 2017 (diana.c.dima@gmail.com)    

classdef lcmv_args
    
    properties
        
        normalize = true;
        marker = 'onset';
        trialfun = 'ft_trialfun_general';
        lowpass = 100;
        bandpass = [];
        prestim = 0.5;
        poststim = 0.7;
        baseline = [-0.5 0];
        resamplefs = 300;
        fixedori = true;
        plot = 0;
        mnn = false;
        sourcemodel = 'standard_sourcemodel3d10mm';
        
    end
    
    methods
        
        function obj = lcmv_args(varargin)
            
            if ~isempty(varargin)
                obj.normalize = varargin{1};
                obj.marker = varargin{2};
                obj.lowpass = varargin{3};
                obj.bandpass = varargin{4};
                obj.prestim = varargin{5};
                obj.poststim = varargin{6};
                obj.baseline = varargin{7};
                obj.resamplefs = varargin{8};
                obj.fixedori = varargin{9};
                obj.plot = varargin{10};
                obj.mnn = varargin{11};
                obj.sourcemodel = varargin{12};
            end;
            
        end
        
    end
    
end
        
