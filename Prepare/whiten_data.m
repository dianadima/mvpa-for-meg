function [mnn_train_data,mnn_test_data] = whiten_data(train_data,train_labels,test_data)
% multivariate noise normalization for binary MEG data (see Guggenmos et al. 2018, NeuroImage)
% data is sensor x time x trial for 2 conditions, divided into training and test sets
%
% DC Dima 2018 (diana.c.dima@gmail.com)

if size(train_data,1)~=size(test_data,1) || size(train_data,2)~=size(test_data,2)
    error('Training and test data must have the same number of features (channels/sources/timepoints).')
end

class_id = unique(train_labels); 
class1 = train_labels==class_id(1); 
class2 = train_labels==class_id(2); %this needs to be done separately for each condition

sigma_time = zeros(2,size(train_data,2), size(train_data,1), size(train_data,1));
for t = 1:size(train_data,2)
    sigma_time(1,t,:,:) = cov1para(squeeze(train_data(:,t,class1))');
    sigma_time(2,t,:,:) = cov1para(squeeze(train_data(:,t,class2))');
end

sigma_time = squeeze(mean(sigma_time,1)); %average across conditions
sigma_inv = (squeeze(mean(sigma_time,1)))^-0.5;

mnn_train_data = nan(size(train_data));
mnn_test_data = nan(size(test_data));

for t = 1:size(train_data,2)
    mnn_train_data(:,t,:) = (squeeze(train_data(:,t,:))'*sigma_inv)';
    mnn_test_data(:,t,:) = (squeeze(test_data(:,t,:))'*sigma_inv)';
end



end

