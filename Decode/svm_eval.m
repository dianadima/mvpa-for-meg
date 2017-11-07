classdef svm_eval
    
    properties (Access = public)
        
        %here are the defaults
        AUC = false;
        ind_folds = false;
        weights = false;
        plotROC = false;
        svm_model = false;
       
    end
    
    methods
        
        function obj = svm_eval(varargin)
            
            if ~isempty(varargin)
                
                obj.AUC = varargin{1};
                obj.plotROC = varargin{2};
                obj.weights = varargin{3};
                obj.ind_folds = varargin{4};
                obj.svm_model = varargin{5};
            end;
        end
                
    end
    
end