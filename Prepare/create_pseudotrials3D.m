function [ newdata, newlabels ] = create_pseudotrials3D( data, labels, n_trials, n_perm)
% Creates averages of trials to increase SNR for multivariate analysis.
% (!) Make sure that last dimension of data is trial dimension. No singleton dimensions.
% Inputs:
%       data = 3-dimensional matrix with trials as last dimension.
%       labels = class labels, of same length as last dimension of data. Any number of classes is supported. For one class, input labels as [].
%       n_trials = number of trials to average together (e.g., 5).
%       n_perm = number of times to repeat the averaging with random assignment of trials to subgroups.
%                 If n_perm > 1, output will contain permutations x data x trials.
% Regardless of the number of permutations, it returns data (features in original n dimensions x pseudotrials) and labels (pseudotrials x 1 vector). 
%
% DC Dima 2018 (diana.c.dima@gmail.com)

nchn = size(data,1);
nwin = size(data,2);

classes = unique(labels);
ncon = length(classes);
permidx = cell(ncon,1);
condidx = cell(ncon,1);
for i = 1:ncon
    idx = find(labels==classes(i));
    condidx{i} = idx;
    
    pidx = nan(numel(idx), n_perm);
    for p = 1:n_perm
        nidx = 1:numel(idx);
        pidx(:,p) = nidx(randperm(numel(idx))); 
    end
    permidx{i} = pidx;
end

newlabels = [];
newdata = [];

for c = 1:ncon
    
        tmpdata = data(:,:,condidx{c}); 
        pidx = permidx{c}; %nperm x trl
        
        for p = 1:n_perm
            
            tmpd = tmpdata(:,:,pidx(:,p));
            tidx = 1:n_trials:size(tmpd,3)-n_trials;
            tmp = nan(size(tmpd,1),size(tmpd,2),numel(tidx));
            for t = 1:numel(tidx)
                tmp(:,:,t) = nanmean(tmpd(:,:,tidx(t):tidx(t)+n_trials-1),3);
            end
            newdata = cat(3, newdata, tmp);
            newlabels = [newlabels; classes(c)*ones(size(tmp,3),1)];
        end
        
    
end

end