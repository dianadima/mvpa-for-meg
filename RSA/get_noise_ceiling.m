function [ noise_ceiling ] = get_noise_ceiling( subject_rdms)
% Calculate upper and lower bounds of noise ceiling using leave-one-out approach (Nili et al., 2014)
% Input: subject rdm matrix (features x subjects)
% Output: noise ceiling bounds (upper,lower)
%
% DC Dima 2018 (diana.c.dima@gmail.com)

average_rdm = squeeze(nanmean(subject_rdms,2)); %average last dim
nc_sub = zeros(2,size(subject_rdms,2));


for idx = 1:size(subject_rdms,2)
    
    sub_rdm = squeeze(subject_rdms(:,idx));
    ave_lso = subject_rdms; ave_lso(:,idx) = [];
    ave_lso = squeeze(nanmean(ave_lso,2));
    nc_sub(1,idx) = corr(sub_rdm(:), average_rdm(:), 'type', 'Spearman');
    nc_sub(2,idx) = corr(sub_rdm(:), ave_lso(:), 'type', 'Spearman');
    
end;

noise_ceiling = nanmean(nc_sub,2);   

end

