function [ correlation_matrix ] = time_resolved_rsa (meg_rdm, models)
%UNTITLED6 Summary of this function goes here
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
    meg_rdm = reshape(meg_rdm, size(meg_rdm,1), size(meg_rdm,2), 1,1);
elseif length(size(meg_rdm))==3
    meg_rdm = reshape(meg_rdm, size(meg_rdm,1), size(meg_rdm,2), size(meg_rdm,3),1);
end;

%prepare matrix of models ready for correlation
mmat = zeros(size(models,1)*size(models,2),size(models,3));
for m = 1:size(models,3)
    
    model = squeeze(models(:,:,m));
    mmat(:,m) = model(:);
    
end;

%here we correlate stuff...
correlation_matrix = zeros(size(meg_rdm,3), size(meg_rdm,4), size(models,3));

for c = 1:size(meg_rdm,4) %in space
    for t = 1:size(meg_rdm,3) %in time
        
        rdm = squeeze(meg_rdm(:,:,t,c));
        x = [rdm(:) mmat];
        
        corrmat = partialcorr(x, 'type', 'Spearman');
        correlation_matrix(t,c,:) = corrmat(1,2:size(models,3)+1);
        
    end;
    
end;

correlation_matrix = squeeze(correlation_matrix);   

end

