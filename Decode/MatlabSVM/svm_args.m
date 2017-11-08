classdef svm_args
    
    properties (Access = public)
        
        %here are the defaults
        kernel = 'linear';
        boxconstraint = 1;
        kfold = 5;
        standardize = true;
        AUC = false;
        ind_folds = false;
        weights = false;
        plotROC = false;
        svm_model = false;
       
    end
    
    methods
        
        function obj = svm_args(varargin)
            
            if ~isempty(varargin)
                obj.kernel = char(varargin{1});
                obj.boxconstraint = varargin{2};
                obj.kfold = varargin{3};
                obj.standardize = varargin{4};
                obj.AUC = varargin{5};
                obj.plotROC = varargin{6};
                obj.weights = varargin{7};
                obj.ind_folds = varargin{8};
                obj.svm_model = varargin{9};                
            end;
        end
                
    end
    
end