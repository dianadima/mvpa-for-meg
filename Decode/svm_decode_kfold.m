function results = svm_decode_kfold (data, labels, varargin)
% Basic function using LibLinear implementation of svm for classification. Only implements kfold crossvalidation.
% Inputs: data(trials x features), labels. 
% Optional: svm parameters and desired optional performance metrics.
%
% Name-value pairs: 
% 'solver', default 1: L2-regularized dual problem solver (LibLinear);
% 'boxconstraint', default 1
% 'kfold', default 5
% 'iterate_cv', default 1, number of cross-validation iterations
% 'standardize', default true (recommended; across training set)
% 'weights', default false (output vector of classifier weights & activation patterns associated with them cf. Haufe 2014)
%            Note: weights are computed by retraining classifier on whole dataset
% 'cv_indices', prespecified training set cross-val indices (cvpartition object or indices x folds matrix/logical, with 1 = training set)
%
% Outputs results structure with non-optional metrics: accuracy, Fscore, sensitivity, specificity. 
% Optional (as above): weights and weight-derived patterns. (Set 'weights',true)
%
% DC Dima 2018 (diana.c.dima@gmail.com)

%parse inputs
svm_par = args.svm_args;
list = fieldnames(svm_par);
p = inputParser;
for ii = 1:length(properties(args.svm_args))
    addParameter(p, list{ii}, svm_par.(list{ii}));
end
parse(p, varargin{:});
svm_par = p.Results;
clear p

if abs(nargin)<2 
    error('Data and labels are needed as inputs.')
end

if ~isa(data, 'double')
    data = double(data);
end

labels = labels(:);
results = struct;

for icv = 1: svm_par.iterate_cv

    %cross-validation indices - convoluted but can be kept constant if need be
    if isempty(svm_par.cv_indices)
        cv = cvpartition(labels, 'kfold', svm_par.kfold); %stratified by default
    else
        if iscell(svm_par.cv_indices)
            cv = svm_par.cv_indices{icv};
        else
            cv = svm_par.cv_indices;
        end
    end
    
    if isa(cv,'cvpartition')
        cv_train = zeros(length(labels),1); cv_test = cv_train;
        for ii = 1:svm_par.kfold
            cv_train(:,ii) = cv.training(ii);
            cv_test(:,ii) = cv.test(ii);
        end
    else
        if size(cv,2)~=svm_par.kfold
            error('Crossval indices must be supplied in indices x folds matrix or logical array.');
        end
        cv_train = cv;
        cv_test = abs(cv_train-1);
    end
    
    cv_train = logical(cv_train); cv_test = logical(cv_test);
    cv_idx = cv_train; %this will be saved for later
    allscore = zeros(length(labels),1); accuracy = zeros(3,svm_par.kfold);
    
    for ii = 1:svm_par.kfold
        
        %scale values using range and minimum of training set
        kdata = data; %ensures we keep original data
        if svm_par.standardize
            data = (kdata - repmat(min(data(cv_train(:,ii),:), [], 1), size(kdata, 1), 1)) ./ repmat(max(data(cv_train(:,ii),:), [], 1) - min(data(cv_train(:,ii),:), [], 1), size(kdata, 1), 1);
        end
        data(isnan(data)|isinf(data)) = 0;
        svm_model = train(labels(cv_train(:,ii)), sparse(data(cv_train(:,ii),:)), sprintf('-s %d -c %d -q 1', svm_par.solver, svm_par.boxconstraint)); %dual-problem L2 solver with C=1
        [allscore(cv_test(:,ii)), accuracy(:,ii), ~] = predict(labels(cv_test(:,ii)), sparse(data(cv_test(:,ii),:)), svm_model, '-q 1');
        
    end

    results.Accuracy(icv) = mean(accuracy(1,:));
    results.AccuracyMSError(icv) = mean(accuracy(2,:));
    results.AccuracyFold(icv,:) = accuracy(1,:);
    results.Confusion{icv} = confusionmat(labels,allscore);
    results.Sensitivity(icv) = results.Confusion{icv}(1,1)/(sum(results.Confusion{icv}(1,:))); %TP/allP = TP/(TP+FN)
    results.Specificity(icv) = results.Confusion{icv}(2,2)/(sum(results.Confusion{icv}(2,:))); %TN/allN = TN/(FP+TN)
    PP = results.Confusion{icv}(1,1)/(sum(results.Confusion{icv}(:,1))); %positive predictive value: class1
    NP = results.Confusion{icv}(2,2)/(sum(results.Confusion{icv}(:,2))); %negative predictive value: class2
    results.Fscore1(icv) = (2*PP*results.Sensitivity(icv))/(PP+results.Sensitivity(icv));
    results.Fscore2(icv) = (2*NP*results.Specificity(icv))/(NP+results.Specificity(icv));
    results.WeightedFscore(icv) = ((sum(results.Confusion{icv}(:,1))/sum(results.Confusion{icv}(:)))*results.Fscore1(icv)) + ((sum(results.Confusion{icv}(:,2))/sum(results.Confusion{icv}(:)))*results.Fscore2(icv));
    results.cv_indices(icv,:,:) = cv_idx; %this can be reused
    results.PredictedLabels(:,icv) = allscore;
end

%calculate weights and compute activation patterns as per Haufe (2014)
if svm_par.weights
    data = (kdata - repmat(min(kdata, [], 1), size(kdata, 1), 1)) ./ repmat(max(kdata, [], 1) - min(kdata, [], 1), size(kdata, 1), 1);
    svm_model = train(labels, sparse(data), sprintf('-s %d -c %d -q 1', svm_par.solver, svm_par.boxconstraint));
    results.Weights = svm_model.w;
    results.WeightPatterns = abs(cov(data)*results.Weights'/cov(data*results.Weights'));
    results.WeightPatternsNorm = (results.WeightPatterns-min(results.WeightPatterns))/(max(results.WeightPatterns)-min(results.WeightPatterns));
end


end
        
