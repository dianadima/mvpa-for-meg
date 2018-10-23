function [ newdata, newlabels ] = create_pseudotrials( data, labels, n_trials, n_perm)
% Creates averages of trials to increase SNR for multivariate analysis.
% (!) Make sure that last dimension of data is trial dimension.
% Inputs:
%       data = n-dimensional matrix with trials as last dimension.
%       labels = class labels, of same length as last dimension of data. Any number of classes is supported.
%       n_trials = number of trials to average together (e.g., 5).
%       n_perm = number of times to repeat the averaging with random assignment of trials to subgroups.
%                 If n_perm > 1, output will contain permutations x data x trials.
% If number of permutations is >1, it returns data (features x trials x
% permutations) and labels (trials x permutations). Otherwise returns data
% (features x trials) and labels (trials column vector).
%
trldim = ndims(data); %trials is last dimension, usually convenient for fieldtrip
sz = [trldim 1:trldim-1]; data = permute(data, sz); %move trials to first dimension
sz = size(data); %store new size

if isempty(labels)
    n_cond = 1;
    idx = repmat({':'}, 1, ndims(data));
    idx{1,1} = 1:size(data, 1);
    trlsize = floor(size(data,1)/n_trials);
else
    n_cond = length(unique(labels));
    classes = unique(labels);
    idx = repmat({':'}, n_cond, ndims(data));
    trlsize = zeros(1,n_cond);
    for i = 1:n_cond
        idx{i,1} = find(labels==classes(i));
        if length(idx{i,1})>1
            trlsize(i) = floor(length(find(labels==classes(i)))/n_trials);
        else
            trlsize(i) = length(idx{i,1});
        end;
    end;
    trlsize = sum(trlsize);
end;

newsz = [n_perm trlsize sz(2:end)];
newdata = zeros(newsz);
newlabels = cell(1,n_perm);
    
for p = 1:n_perm
    
    condata = [];
    
    for c = 1:n_cond
    
        tmpdata = data(idx{c,:});
        
        if size(tmpdata,1)>1 %do we have more than 1 trial?
            
            tmpdata = shuffle(tmpdata, 1);
            tp = 1:n_trials:size(tmpdata, 1)-n_trials+1;
            tmpsz = [length(tp) sz(2:end)];
            avedata = zeros(tmpsz);
            tmpidx = idx;
            
            for t = 1:length(tp)
                
                tmpidx{c,1} = tp(t):tp(t)+n_trials-1;
                pseudo = tmpdata(tmpidx{c,:});
                pseudo = squeeze(nanmean(pseudo,1));
                avedata(t,:) = pseudo(:);
                
            end;
            
        else
            
            avedata = tmpdata;
            
        end;
        
        condata = cat(1, condata, avedata);
                
        newlabels{p} = [newlabels{p}; c*ones(size(avedata,1),1)]; 
    
    end;
    
    newdata(p,:,idx{1,2:end}) = condata;
    
end;

%reshape into data*pseudotrials (over permutations)
newsz = 1:length(newsz); newsz(length(newsz)+1) = 2; newsz(length(newsz)+1) = 1; newsz(1:2)=[]; 
newdata = permute(newdata, newsz); %trials are last dim
%newdata = squeeze(newdata); %one permutation case
newlabels = cat(2,newlabels{:});


end