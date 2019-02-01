function [ dynamic_noise_ceiling ] = time_resolved_noise_ceiling( subject_rdms )
%Calculate time and space - resolved upper and lower bounds of noise ceiling using leave-one-out approach (Nili et al., 2014)
%Input: subject rdm matrix (features x time/space x subjects, or features x time x space x subjects)
%Output: noise ceiling bounds (upper,lower)
%
%DC Dima 2018 (diana.c.dima@gmail.com)

%deal with one time window case
if length(size(subject_rdms))==3
    subject_rdms = reshape(subject_rdms, size(subject_rdms,1), size(subject_rdms,2), 1,1,size(subject_rdms,3));
elseif length(size(subject_rdms))==4
    subject_rdms = reshape(subject_rdms, size(subject_rdms,1), size(subject_rdms,2), size(subject_rdms,3),1,size(subject_rdms,4));
end;

dynamic_noise_ceiling = zeros(size(subject_rdms,3), size(subject_rdms,4), 2);

for c = 1:size(subject_rdms,4) %in space
    for t = 1:size(subject_rdms,3) %in time
        
        rdm = squeeze(subject_rdms(:,:,t,c,:));
        
        dynamic_noise_ceiling(c,t,:) = get_noise_ceiling(rdm);

        
    end;
    
end;

dynamic_noise_ceiling = squeeze(dynamic_noise_ceiling);   

end