classdef svm_props
    
    properties
        
        args = svm_args;
        eval = svm_eval;
        
    end
    
    methods
        
        function obj = svm_props(args,eval)
            
            obj.args = args;
            obj.eval = eval;
            
        end
        
    end
    
end