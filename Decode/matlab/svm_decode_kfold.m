function results = svm_decode_kfold (data, labels, varargin)
% Basic function using MATLAB implementation of svm for classification. Only implements kfold crossvalidation.
% Inputs: data(trials x features), labels. Optional: svm
% parameters and desired optional performance metrics.
% Name-value pairs: 
% 'kernel', default 'linear'
% 'boxconstraint', default 1
% 'kfold', default 5
% 'standardize', default true (recommended; across training & test sets)
% 'AUC', default false (return area under curve; note - can increase computation time, as need to convert binary SVM scores to posterior probabilities)
% 'ind_folds', default false (return accuracy per fold as well as the averaged accuracy across folds)
% 'weights', default false (output vector of classifier weights & activation patterns associated with them cf. Haufe 2014)
% 'plotROC', false (plot receiver operating characteristic curve)
% 'svm_model', false (keep trained svm model)
%
% Outputs results structure with non-optional metrics: accuracy, Fscore, sensitivity, specificity. Structure can be accessed using e.g. accuracy = cell2mat({results.Accuracy}).
% Optional (as above): per-fold accuracy, weights and weight-derived patterns, SVM model, AUC and plot of AUROC.
%
% DC Dima 2017 (diana.c.dima@gmail.com)

%parse inputs
svm_par = svm_args_matlab;
list = fieldnames(svm_par);
p = inputParser;
for i = 1:length(properties(svm_args_matlab))
    addParameter(p, list{i}, svm_par.(list{i}));
end;
parse(p, varargin{:});
svm_par = p.Results;
clear p;

if abs(nargin)<2 
    error('Data and labels are needed as inputs.')
end;

%classification and main results
cp = classperf (labels);
SVM_model = fitcsvm(data, labels, 'KernelFunction', svm_par.kernel, 'BoxConstraint', svm_par.boxconstraint, 'Crossval', 'on', 'KFold', svm_par.kfold, 'Standardize', svm_par.standardize); %train
[pred,scores] = kfoldPredict(SVM_model); %test
cp = classperf(cp, pred); %update performance-tracking object

%non-optional outputs
results.Accuracy = cp.CorrectRate;
results.Fscore  = ((2*((cp.PositivePredictiveValue*cp.Sensitivity)/(cp.PositivePredictiveValue + cp.Sensitivity))) + (2*((cp.NegativePredictiveValue*cp.Specificity)/(cp.NegativePredictiveValue + cp.Specificity))))/2;
results.Confusion = cp.CountingMatrix;
results.Sensitivity = cp.Sensitivity;
results.Specificity = cp.Specificity;

%optional outputs
if svm_par.AUC
    svm_post = fitSVMPosterior(SVM_model);
    [Xsvm,Ysvm,~,results.AUC] = perfcurve(labels, scores(:,svm_post.ClassNames(1)),'1');
end;

if svm_par.plotROC
    figure;
    plot(Xsvm(:,1),Ysvm(:,1), 'k');
    xlabel('False positive rate');
    ylabel('True positive rate');
    title('ROC Curve', 'FontWeight', 'normal');
end;

if svm_par.weights
    final_SVM = fitcsvm(data, labels, 'KernelFunction', svm_par.kernel, 'BoxConstraint', svm_par.boxconstraint, ...
        'Standardize', svm_par.standardize, 'IterationLimit', 1e8, 'CacheSize', 'maximal');
    %get feature weights (Betas/coeffs) - remember to use absolute values or squares
    results.Weights = final_SVM.Beta';
    % get patterns as per Haufe 2014
    results.WeightPatterns = abs(results.Weights*cov(data));
end;

if svm_par.ind_folds
    results.FoldAccuracy = 1 - kfoldLoss(svm_post,'mode','individual');
end;

if svm_par.svm_model
    results.SVMmodel = compact(SVM_model);
    results.SVMpost = compact(svm_post);
end;

end
        
