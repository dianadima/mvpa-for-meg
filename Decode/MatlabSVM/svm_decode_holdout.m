function results = svm_decode_holdout (train_data, train_labels, test_data, test_labels, varargin)
% Inputs: training and testing data(trials x features), training and testing abels. 
% Optional: svm parameters and desired optional performance metrics.
% Name-value pairs: 
% 'kernel', default 'linear'
% 'boxconstraint', default 1
% 'standardize', default true (recommended; across training & test sets)
% 'AUC', default false (return area under curve; note - can increase computation time, as need to convert binary SVM scores to posterior probabilities)
% 'ind_folds', default false (return accuracy per fold as well as the averaged accuracy across folds)
% 'weights', default false (output vector of classifier weights & activation patterns associated with them cf. Haufe 2014)
% 'plotROC', false (plot receiver operating characteristic curve)
% 'svm_model', false (keep trained svm model)
%
% Outputs results structure with non-optional metrics: accuracy, Fscore, sensitivity, specificity. Structure can be accessed using e.g. accuracy = cell2mat({results.Accuracy}).
% Optional (as above): per-fold accuracy, weights and weight-derived patterns, SVM model, AUC and plot of AUROC.
% Basic function using MATLAB implementation of svm for classification, using independent training and test sets.

%parse inputs
svm_par = svm_args;
list = fieldnames(svm_par);
p = inputParser;
for i = 1:length(properties(svm_args))
    addParameter(p, list{i}, svm_par.(list{i}));
end;
parse(p, varargin{:});
svm_par = p.Results;
clear p;

if abs(nargin)<2 
    error('Data and labels are needed as inputs.')
end;

%classification and main results
cp = classperf (test_labels);
SVM_model = fitcsvm(train_data, train_labels, 'KernelFunction', svm_par.kernel, 'BoxConstraint', 1,  'Standardize', svm_par.standardize); %train
[pred,scores] = predict(SVM_model, test_data); %test
cp = classperf(cp, pred); %update performance-tracking object

%non-optional outputs
results.Accuracy = cp.CorrectRate;
results.Fscores  = ((2*((cp.PositivePredictiveValue*cp.Sensitivity)/(cp.PositivePredictiveValue + cp.Sensitivity))) + (2*((cp.NegativePredictiveValue*cp.Specificity)/(cp.NegativePredictiveValue + cp.Specificity))))/2;
results.Confusion = cp.CountingMatrix;
results.Sensitivity = cp.Sensitivity;
results.Specificity = cp.Specificity;

%optional outputs
if svm_par.AUC==true
    svm_post = fitSVMPosterior(SVM_model);
    [Xsvm,Ysvm,~,results.AUC] = perfcurve(labels, scores(:,svm_post.ClassNames(1)),'1');
end;

if svm_par.plotROC==true
    figure;
    plot(Xsvm(:,1),Ysvm(:,1), 'k');
    xlabel('False positive rate');
    ylabel('True positive rate');
    title('ROC Curve', 'FontWeight', 'normal');
end;

if svm_par.weights==true
    %get feature weights (Betas/coeffs) - remember to use absolute values or squares
    results.Weights = SVM_model.Beta';
    % get patterns as per Haufe 2014
    results.WeightPatterns = abs(results.Weights*cov(train_data));
end;

if svm_par.ind_folds==true
    results.FoldAccuracy = 1 - kfoldLoss(svm_post,'mode','individual');
end;

if svm_par.svm_model==true
    results.SVMmodel = compact(SVM_model);
    results.SVMpost = compact(svm_post);
end;

end
        