function [ accuracy, random_results ] = randomize_svm( data, labels, num_iterations)
% Inputs: data, labels, number of randomization iterations desired (i.e. 1000)
% Outputs: randomized accuracy and (optional) entire results structure for each iteration.
% Classification on randomized data (shuffles labels across training and
% test sets).

accuracy = zeros(1,num_iterations);
random_results = cell(1,num_iterations);

for idx = 1:num_iterations
    
    rand_labels = labels(randperm(length(labels)));
    results = svm_decode(data, rand_labels);
    accuracy(idx) = results.Accuracy;
    random_results{idx} = results;
       
end

