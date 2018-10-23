function [ rand_stat, pvalue, maxvalue ] = randomize_accuracy(accuracy, num_iterations, varargin)
% Sign test on accuracies. Core (scalar accuracy x number of
% observations/subjects).

addParameter(p, 'chance_level', 50);
addParameter(p, 'statistic', 'mean');
parse(p, varargin{:});

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
        pvalue = (length(find(rand_stat>=mean(acc_dm)))+1)/(num_iterations+1);
        maxvalue = max(rand_stat);
        
    case 'tstat'
        
        rand_stat = mean(rand_accuracy,2)/(std(rand_accuracy,[],2)/sqrt(size(accuracy,2)));
        obs_stat = mean(acc_dm)/(std(acc_dm)/sqrt(length(acc_dm)));
        pvalue = (length(find(rand_stat>=obs_stat))+1)/(num_iterations+1);
        maxvalue = max(rand_stat);
        
        
end
     
end

