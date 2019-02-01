function [ random_corr ] = randomize_rsa( meg_rdm, models, num_iterations )
% Run Spearman's correlation between shuffled MEG RDM and models to estimate null empirical distribution.
% The RDMs should be vectores obtained from symmetric matrices using tril or triu (lower or upper triangles WITHOUT the diagonal).
%
% Inputs: meg_rdm (features x time), models (features x time x models)
%         OR meg_rdm (features x space x time), models (features x space x time x models)
%
% Output: randomized correlation coefficients (iterations x space/time x models)
%
% DC Dima 2017 (diana.c.dima@gmail.com)

if ismatrix(meg_rdm)
    
    random_corr = zeros(num_iterations, size(meg_rdm,2), size(models,3));
    
    for idx = 1:num_iterations
        
        meg_rdm = meg_rdm(randperm(size(meg_rdm,1)),:); %shuffle RDM - same scheme across time
        
        for t = 1:size(meg_rdm,2)
            
            rdm = squeeze(meg_rdm(:,t));
            m = squeeze(models(:,t,:));
            x = [rdm(:) m];
            corrmat = partialcorr(x, 'type', 'Spearman');
            random_corr(idx,t,:) = corrmat(1,2:size(models,3)+1);
            
        end
        
    end
    
elseif ndims(meg_rdm)==3
    
    random_corr = zeros(num_iterations, size(meg_rdm,2), size(meg_rdm,3), size(models,4));
    
    for idx = 1:num_iterations
        
        meg_rdm = meg_rdm(randperm(size(meg_rdm,1)),:,:); %shuffle RDM - same scheme across time & space
        
        for s = 1:size(meg_rdm,2)
            
            for t = 1:size(meg_rdm,3)
                
                rdm = squeeze(meg_rdm(:,s,t));
                m = squeeze(models(:,s,t,:));
                x = [rdm(:) m];                
                corrmat = partialcorr(x, 'type', 'Spearman');
                random_corr(idx,s,t,:) = corrmat(1,2:size(models,3)+1);
                
            end
        end
        
    end
    
end

