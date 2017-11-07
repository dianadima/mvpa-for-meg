classdef svm_args
    
    properties (Access = public)
        
        %here are the defaults
        kernel = 'linear';
        boxconstraint = 1;
        kfold = 5;
        standardize = true;
       
    end
    
    methods
        
        function obj = svm_args(varargin)
            
            if ~isempty(varargin)
                obj.kernel = char(varargin{1});
                obj.boxconstraint = varargin{2};
                obj.kfold = varargin{3};
                obj.standardize = varargin{4};
            end;
        end
                
    end
    
end