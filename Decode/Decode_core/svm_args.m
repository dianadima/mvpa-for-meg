classdef svm_args
    
    properties (Access = public)
        
        %here are the defaults
        solver = 1; %liblinear: 1: L2 dual-problem; 2: L2 primal; 3:L2RL1L...
        boxconstraint = 1;
        kfold = 5;
        cv_indices = [];
        iterate_cv = 1;
        standardize = true;
        mnn = true;
        weights = false;
        kernel = 0; %libSVM: 	0: linear; 1: polynomial; 2:radial basis function; 3:sigmoid
       
    end
    
    methods
        
        function obj = svm_args(varargin)
            
            if ~isempty(varargin)
                obj.solver = varargin{1};
                obj.boxconstraint = varargin{2};
                obj.kfold = varargin{3};
                obj.cv_indices = varargin{4};
                obj.iterate_cv = varargin{5};
                obj.standardize = varargin{6};
                obj.mnn = varargin{7};
                obj.weights = varargin{8};
                obj.kernel = varargin{9};
            end;
        end
                
    end
    
end