function [ accuracy, random_results ] = randomize_svm( data, labels, num_iterations, varargin)
% Inputs: data, labels, number of randomization iterations desired (i.e. 1000) + svm decoding parameters (for kfold crossvalidation).
% If you want to randomize holdout accuracy, the last 2 arguments should be a test set and test labels.
% Outputs: randomized accuracy and (optional) entire results structure for each iteration.
% Classification on randomized data (shuffles labels across training and test sets).

accuracy = zeros(1,num_iterations);
random_results = cell(1,num_iterations);

for idx = 1:num_iterations
    
    rand_labels = labels(randperm(length(labels)));
    
    if isempty(varargin)
        results = svm_decode_kfold(data, rand_labels);
    elseif length(varargin)==2
        test_labels = varargin{2}; test_labels = test_labels(randperm(length(test_labels)));
        results = svm_decode_holdout(data, rand_labels, varargin{1}, test_labels);
    else
        error('Wrong number of arguments. Read help file');
    end;
    
    accuracy(idx) = results.Accuracy;
    random_results{idx} = results;
       
end

