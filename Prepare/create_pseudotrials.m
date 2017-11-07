function [ newdata ] = create_pseudotrials( data, n_trials, n_perm, varargin )
%Creates averages of trials to increase SNR for multivariate analysis.
% Make sure that last dimension of data is trial dimension.
% Optional parameter: cell array containing indices of trials
% belonging to each condition. Example {1:40, 41:80}.

if ~isempty(varargin)
    conditions = varargin{1};
    n_cond = length(varargin{1});    
else
    n_cond = 1;
end;

trldim = ndims(data); %trials is last dimension
newsz = [n_cond n_perm floor(size(data,trldim)/n_trials) size(data)];
newdata = zeros(newsz);
idx = repmat({':'},1,ndims(data)); %how we index data to be averaged

for c = 1:n_cond
    
    for p = 1:n_perm
        
        if n_cond>1
            
            idx{trldim} = conditions{c};
            tmpdata = data(idx{:});
            tmpdata = shuffle(tmpdata, trldim);
            tp = 1:n_trials:length(conditions{c})-n_trials+1;
            
        else
            tmpdata = shuffle(data,trldim);
            tp = 1:n_trials:size(data,trldim)-n_trials+1;
            
        end;
        
        for t = 1:length(tp)
            
            idx{trldim} = tp(t):tp(t)+n_trials-1;
            pseudo = tmpdata(idx{:});
            pseudo = squeeze(nanmean(pseudo,trldim));
            newdata(c,p,t,:) = pseudo(:);
            
        end;
    end;
    
end;

%reshape into data*pseudotrials (over permutations)
newsz = size(newdata); newsz(2) = newsz(2)*newsz(3); newsz(3) = [];
newdata = reshape(newdata,newsz);
newsz = 1:length(newsz); newsz(length(newsz)+1) = 2; newsz(2)=[]; 
newdata = permute(newdata, newsz); %trials are last dim
newdata = squeeze(newdata); %one condition case


end