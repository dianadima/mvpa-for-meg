function [ rand_stat, obs_stat, pvalue ] = randomize_accuracy(accuracy, varargin)
% Sign test on accuracies. 
% Inputs: accuracy (length is number of observations/subjects).
% Name-value optional inputs: 'num_iterations' (default 5000) - number of randomizations;
%                             'chance_level' (default 50) - will be subtracted before sign flipping;
%                             'statistic' (default 'mean') - which statistic to randomize, 'mean' or 'tstat'
% Outputs: rand_stat contains randomized statistic (length is num_iterations);
%          p-value (one-tailed: number of randomizations exceeding observed statistic).
%
% DC Dima 2018 (diana.c.dima@gmail.com)

p = inputParser;
addParameter(p, 'num_iterations',5000);
addParameter(p, 'chance_level', 50);
addParameter(p, 'statistic', 'mean');
parse(p, varargin{:});

num_iterations = p.Results.num_iterations;
acc_dm = accuracy-p.Results.chance_level; %demean accuracy
rand_sign = sign(randn(num_iterations,length(acc_dm)));
if size(acc_dm,1) == 1
    rand_accuracy = repmat(acc_dm, num_iterations,1).*rand_sign;
elseif size(acc_dm,2) == 1
    rand_accuracy = repmat(acc_dm, 1, num_iterations)'.*rand_sign; %iterations x subjects
elseif ~isvector(acc_dm)
    error('Only vector data allowed')
end
    
switch p.Results.statistic
    
    case 'mean'
        
        rand_stat = mean(rand_accuracy,2);
        obs_stat = mean(acc_dm);
        pvalue = (length(find(rand_stat>=mean(acc_dm)))+1)/(num_iterations+1);
        
    case 'tstat'
        
        rand_stat = mean(rand_accuracy,2)/(std(rand_accuracy,[],2)/sqrt(size(accuracy,2)));
        obs_stat = mean(acc_dm)/(std(acc_dm)/sqrt(length(acc_dm)));
        pvalue = (length(find(rand_stat>=obs_stat))+1)/(num_iterations+1);       
        
end
     
end

