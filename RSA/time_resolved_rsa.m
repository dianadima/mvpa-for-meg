function [ correlation_matrix ] = time_resolved_rsa (meg_rdm, models)
% Run Representational Similarity RSA by comparing a MEG RDM to models over time.
% Input: meg_rdm (features x time), models (features x time x models)
% Output: correlation matrix (time x models)
%
% DC Dima 2017 (diana.c.dima@gmail.com)

if ismatrix(meg_rdm) %features by time
    
    correlation_matrix = zeros(size(meg_rdm,2), size(models,3));
    
    for t = 1:size(meg_rdm,3) %in time
        
        rdm = squeeze(meg_rdm(:,t));
        m = squeeze(models(:,t,:));
        x = [rdm(:) m];
        corrmat = partialcorr(x, 'type', 'Spearman');
        correlation_matrix(t,:) = corrmat(1,2:size(models,3)+1);
        
    end
    
    
elseif ndims(meg_rdm)==3 %features by space by time
    
    correlation_matrix = zeros(size(meg_rdm,2), size(meg_rdm,3), size(models,4));
    
    for s = 1:size(meg_rdm,2) %in space
        
        for t = 1:size(meg_rdm,3) %in time
            
            rdm = squeeze(meg_rdm(:,s,t));
            m = squeeze(models(:,s,t,:));
            x = [rdm(:) m];
            corrmat = partialcorr(x, 'type', 'Spearman');
            correlation_matrix(s,t,:) = corrmat(1,2:size(models,3)+1);
            
        end
        
    end
    
end

correlation_matrix = squeeze(correlation_matrix);

end

