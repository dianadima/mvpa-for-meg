function [ random_corr ] = randomize_rsa( meg_rdm, models, num_iterations )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
if isvector(meg_rdm) || isvector(models)
    error('RDMs need to be supplied in matrix format');
end;

%deal with one model case
if length(size(models))==2
    models = reshape(models, size(models,1), size(models,2), 1);
end;

%deal with one time window case
if length(size(meg_rdm))==2
    meg_rdm = reshape(meg_rdm, size(meg_rdm,1), size(meg_rdm,2), 1);
end;

%here we correlate stuff...
random_corr = zeros(num_iterations, size(meg_rdm,3), size(models,3));

for idx = 1:num_iterations
    
    meg_rdm = meg_rdm(randperm(size(meg_rdm,1)),:,:); %shuffle RDM - same scheme across time
    
    for t = 1:size(meg_rdm,3)
        
        rdm = squeeze(meg_rdm(:,:,t));
        x = rdm(:);
        
        for m = 1:size(models,3)
            
            model = squeeze(models(:,:,m));
            x(:,m+1) = model(:);
            
        end;
        
        corrmat = partialcorr(x, 'type', 'Spearman');
        random_corr(idx,t,:) = corrmat(1,2:size(models,3)+1);
        
    end;
    
end;

end

