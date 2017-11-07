function [ noise_ceiling ] = get_noise_ceiling( subject_rdms)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

average_rdm = squeeze(nanmean(subject_rdms,3)); %average last dim
nc_sub = zeros(2,size(subject_rdms,3));


for idx = 1:size(subject_rdms,3)
    
    sub_rdm = squeeze(subject_rdms(:,:,idx));
    ave_lso = subject_rdms; ave_lso(:,:,idx) = [];
    ave_lso = squeeze(nanmean(ave_lso,3));
    nc_sub(1,idx) = corr(sub_rdm(:), average_rdm(:), 'type', 'Spearman');
    nc_sub(2,idx) = corr(sub_rdm(:), ave_lso(:), 'type', 'Spearman');
    
end;

noise_ceiling = nanmean(nc_sub,2);   

end

