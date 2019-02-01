function [ rand_stat, true_stat, pvalue ] = randomize_labels(predicted_labels, true_labels, varargin)
% Comparing predicted labels to a random distribution of shuffled labels.
%
% Inputs: predicted_labels: column vector or matrix (number observations x number randomizations, e.g. from iterated cross-validation). 
%                           In latter case, each column will be compared to randomized labels to get a mean statistic for each shuffling iteration.
%         true_labels: column vector (labels are scalars, e.g. -1, 1, 2, for 2 classes)
%
% Name-value optional inputs: 'num_iterations' (default 5000) - number of randomizations;
%                             'statistic' (default 'accuracy') - which
%                             statistic to randomize, 'accuracy' or 'sens_spec' (sensitivity & specificity, returned in this order)
% 
% Outputs: rand_stat contains randomized statistic: accuracy x num_iterations, or [sensitivity specificity] x num_iterations
%          true_stat, calculated from the true labels: accuracy or [sensitivity specificity] . 
%                     If several predicted label vectors are given, the mean accuracy/sensitivity & specificity is returned.
%          p-value (one-tailed: number of randomizations exceeding observed statistic).
%
% DC Dima 2019 (diana.c.dima@gmail.com)

p = inputParser;
addParameter(p, 'num_iterations',5000);
addParameter(p, 'statistic', 'accuracy'); % or 'sens_spec'
parse(p, varargin{:});

nperm = p.Results.num_iterations;
ncv = size(predicted_labels,2);
nobs = size(predicted_labels,1);

if length(true_labels)~=nobs
    error('True_labels is a column vector with length equal to number of rows in predicted_labels.')
end

truelabelmat = repmat(true_labels,1,ncv); %where ncv>1, makes a matrix to compare to multiple predicted label vectors

if strcmp(p.Results.statistic, 'accuracy')
    
    rand_stat = nan(1,nperm);
    
    for i = 1:nperm
        
        randlabelmat = truelabelmat(randperm(nobs),:);
        rand_stat(i) = mean(sum(predicted_labels==randlabelmat,1)/nobs);
          
    end

    true_stat = mean(sum(predicted_labels==truelabelmat,1)/nobs);
    pvalue = (length(find(rand_stat>=true_stat))+1)/(nperm+1);
    
elseif strcmp(p.Results.statistic, 'sens_spec')
    
    rand_stat = nan(2,nperm);
    lval = unique(true_labels);
    
    for i = 1:nperm
        
        randlabelmat = truelabelmat(randperm(nobs),:);
        truepos = predicted_labels==randlabelmat&predicted_labels==lval(1);
        trueneg = predicted_labels==randlabelmat&predicted_labels==lval(2);
        rand_stat(1,i) = nanmean(sum(truepos,1)/sum(randlabelmat(:,1)==lval(1)));
        rand_stat(2,i) = nanmean(sum(trueneg,1)/sum(randlabelmat(:,1)==lval(2)));

    end

    truepos = predicted_labels==truelabelmat&predicted_labels==lval(1);
    trueneg = predicted_labels==truelabelmat&predicted_labels==lval(2);
    true_stat(1) = nanmean(sum(truepos,1)/sum(true_labels==lval(1)));
    true_stat(2) = nanmean(sum(trueneg,1)/sum(true_labels==lval(2)));

    pvalue(1) = (length(find(rand_stat(1,:)>=true_stat(1)))+1)/(nperm+1);
    pvalue(2) = (length(find(rand_stat(2,:)>=true_stat(2)))+1)/(nperm+1);
    
else
    error('Statistic not supported, options are accuracy and sens_spec.')
end
        
end