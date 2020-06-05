function [ rand_stat, obs_stat, pvalue ] = randomize_accuracy(accuracy, varargin)
% Sign test on accuracies.
% Inputs: accuracy: vector, 2D or 3D matrix: subjects, subjects x feat, subjects x feat x feat
%                   null distribution estimated consistently across timepoints/ROIs
% Name-value optional inputs: 'num_iterations' (default 5000) - number of randomizations;
%                             'chance_level' (default 50) - will be subtracted before sign flipping;
%                             'statistic' (default 'mean') - which statistic to randomize, 'mean' or 'tstat'
% Outputs: rand_stat contains randomized statistic (num_iterations x statistic);
%          p-value (one-tailed: number of randomizations exceeding observed statistic).
%
% DC Dima 2018 (diana.c.dima@gmail.com)

p = inputParser;
addParameter(p, 'num_iterations',5000);
addParameter(p, 'chance_level', 50);
addParameter(p, 'statistic', 'mean');
parse(p, varargin{:});

num_dim = ndims(accuracy);
num_obs = size(accuracy,1);
num_iterations = p.Results.num_iterations;

acc_dm = accuracy-p.Results.chance_level; %demean accuracy

if isvector(acc_dm), acc_dm = acc_dm(:); end
rand_sign = sign(randn(num_iterations,num_obs));

if num_iterations<=5000 && (num_dim<3 || sum(size(acc_dm))<3000) %vectorize if array not too big
    
    switch num_dim
        case 1
            rand_accuracy = repmat(acc_dm, 1, num_iterations)'.*rand_sign; %iterations x subjects
        case 2
            rand_sign = repmat(rand_sign,[1,1,size(acc_dm,2)]);
            rep_acc = repmat(acc_dm,[1,1,num_iterations]);
            rep_acc = permute(rep_acc,[3 1 2]);
            rand_accuracy = rep_acc.*rand_sign; %it x sub x feat
        case 3
            rand_sign = repmat(rand_sign,[1,1,size(acc_dm,2),size(acc_dm,3)]);
            rep_acc = repmat(acc_dm,[1,1,1,num_iterations]);
            rep_acc = permute(rep_acc,[4 1 2 3]);
            rand_accuracy = rep_acc.*rand_sign; %it x sub x feat x feat
    end
    
    switch p.Results.statistic
        case 'mean'
            rand_stat = squeeze(mean(rand_accuracy,2));
            obs_stat = squeeze(mean(acc_dm,1));
        case 'tstat'
            rand_stat = squeeze(mean(rand_accuracy,2)./(std(rand_accuracy,[],2)./sqrt(num_obs)));
            obs_stat = squeeze(mean(acc_dm,1)./(std(acc_dm,[],1)./sqrt(num_obs)));
    end
    
    
else
    
    sz1 = size(acc_dm,2); sz2 = size(acc_dm,3);
    rand_stat = nan(num_iterations,sz1,sz2);
    switch p.Results.statistic
        case 'mean'
            obs_stat = squeeze(mean(acc_dm,1));
            for i = 1:num_iterations
                rsgn = repmat(squeeze(rand_sign(i,:))',[1 sz1 sz2]);
                racc = acc_dm.*rsgn;
                rand_stat(i,:,:) = mean(racc,1);
            end
        case 'tstat'
            obs_stat = squeeze(mean(acc_dm,1)./(std(acc_dm,[],1)./sqrt(num_obs)));
            for i = 1:num_iterations
                rsgn = repmat(squeeze(rand_sign(i,:)),[1 sz1 sz2]);
                racc = acc_dm.*rsgn;
                rand_stat(i,:,:) = mean(racc,1)./(std(racc,[],1)./(sqrt(num_obs)));
            end
    end
end


if isvector(obs_stat),obs_stat = obs_stat(:); end
pvalue = nan(size(obs_stat));
for i = 1:size(obs_stat,1)
    for ii = 1:size(obs_stat,2)
        pvalue(i,ii) = (length(find(rand_stat(:,i,ii)>obs_stat(i,ii)))+1)/(num_iterations+1);
    end
end



end

