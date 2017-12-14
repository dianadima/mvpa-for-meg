function [ newdata, newlabels ] = create_pseudotrials( data, labels, n_trials, n_perm)
%Creates averages of trials to increase SNR for multivariate analysis.
% Make sure that last dimension of data is trial dimension.
% Optional parameter: cell array containing indices of trials
% belonging to each condition. Example {1:40, 41:80}.

trldim = ndims(data); %trials is last dimension, usually convenient for fieldtrip
sz = [trldim 1:trldim-1]; data = permute(data, sz); %move trials to first dimension
sz = size(data); %store new size

if isempty(labels)
    n_cond = 1;
    idx = repmat({':'}, 1, ndims(data));
    idx{1} = 1:size(data, trldim);
else
    n_cond = length(unique(labels));
    classes = unique(labels);
    idx = repmat({':'}, n_cond, ndims(data));
    for i = 1:n_cond
        idx{i,1} = find(labels==classes(i));
    end;
end;

newsz = [n_perm floor((size(data,1)-1)/n_trials) sz(2:end)];
newdata = zeros(newsz);
newlabels = [];
    
for p = 1:n_perm
    
    condata = [];
    
    for c = 1:n_cond
    
        tmpdata = data(idx{c,:});
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
        
        condata = cat(1, condata, avedata);
                
        if p==1 %no need to do this for each permutation
            newlabels = [newlabels; c*ones(size(avedata,1),1)]; %#ok<AGROW>
        end;
    
    end;
    
    newdata(p,:,idx{1,2:end}) = condata;
    
end;

%reshape into data*pseudotrials (over permutations)
newsz = 1:length(newsz); newsz(length(newsz)+1) = 2; newsz(2)=[]; 
newdata = permute(newdata, newsz); %trials are last dim
newdata = squeeze(newdata); %one permutation case

end