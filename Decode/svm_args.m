classdef svm_args
    
    properties (Access = public)
        
        %here are the defaults
        solver = 1; %liblinear: 1: L2 dual-problem; 2: L2 primal; 3:L2RL1L...
        boxconstraint = 1;
        kfold = 5;
        standardize = true;
        weights = false;
       
    end
    
    methods
        
        function obj = svm_args(varargin)
            
            if ~isempty(varargin)
                obj.solver = varargin{1};
                obj.boxconstraint = varargin{2};
                obj.kfold = varargin{3};
                obj.standardize = varargin{4};
                obj.weights = varargin{5};
            end;
        end
                
    end
    
end