function [train_data,train_labels,test_data,test_labels ] = whiten_data(train_data,train_labels,test_data,test_labels )
%multivariate noise normalization for binary MEG data
%data is sensor x time x trial for 2 conditions

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

for t = 1:size(data,2)
    train_data(:,t,:) = (squeeze(train_data(:,t,:))'*sigma_inv)';
    test_data(:,t,:) = (squeeze(test_data(:,t,:))'*sigma_inv)';
end



end

