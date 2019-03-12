function [ results ] = searchlight_kfold( data, labels, cluster_idx, varargin )
% Performs time-resolved SVM decoding of MEG data, using stratified k-fold cross-validation on each time window, and LibLinear SVM implementation.
% Inputs: data, labels. Data format: channels/sources x time x trials
%         cluster_idx, structure obtained using get_sensor_info or get_source_info - for channel/source selection. You can also just provide numerical indices, in which case you don't need the structure.
% Optional: 
%          'channels', channel set (string or cell array of strings; default: 'MEG').
%          'decoding_window' (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%          'pseudo', default [], create pseudotrials: example [5 100], average groups of 5 trials with 100 random assignments of trials to groups
%          'mnn', default true, perform multivariate noise normalization (recommended, Guggenmos et al. 2018)
%
%           Below are other name-value pairs that control the SVM settings, with defaults:
%           
%          solver = 1; %only applies to liblinear: 1: L2 dual-problem; 2: L2 primal; 3:L2RL1L...
%          boxconstraint = 1; --> C-parameter: Note, we don't have any options for optimizing this, need to write it separately if needed
%          kfold = 5; --> number of folds for k-fold-cross-validation
%          cv_indices = []; --> supply cross-validation indices (e.g. in cvpartition format)
%          standardize = true; --> standardize features using mean and SD of training set (recommended)
%          weights = false; --> calculate weights (by retraining model on whole dataset)
%
% Outputs: structure containing classification performance metrics (for each timepoint and cross-validation round).
%
% DC Dima 2018 (diana.c.dima@gmail.com)

%parse inputs
dec_args = args.decoding_args;
svm_par = args.svm_args;
list = [fieldnames(dec_args); fieldnames(svm_par)];
p = inputParser;

for i = 1:length(properties(args.decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
for ii = i+1:length(properties(args.decoding_args))+length(properties(args.svm_args))
    addParameter(p, list{ii}, svm_par.(list{ii}));
end;
addParameter(p, 'sensor_idx', []);


parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','decoding_window', 'time', 'sensor_idx','pseudo','mnn'}); %converted struct will be fed into decoding function
clear p;

%channel indices for each searchlight iteration
if isstruct(cluster_idx) %the sensor-space case
    chan_idx = arrayfun(@(i) find(ismember({cluster_idx.label},[cluster_idx(i).label; cluster_idx(i).neighblabel])), 1:length(cluster_idx), 'UniformOutput', false); %store all searchlight idx in a cell array
else
    chan_idx = cluster_idx; %the source-space case
end;

%create time axis
if ~isempty(dec_args.time)
    time = dec_args.time;
elseif ~isempty(dec_args.sensor_idx) && isfield(sensor_idx, 'time')
    time = sensor_idx.time;
else
    time = 1:size(data,2);
end;

if length(time)~=size(data,2)
    time = 1:size(data,2);
    fprintf('\nWarning: time axis does not match dataset size. Replacing with default time axis...');
end;

%subselect data corresponding to decoding window requested
if ~isempty(dec_args.decoding_window)
    if ~isempty(find(round(time,3)==dec_args.decoding_window(1),1))
        lims(1) = find(round(time,3)==dec_args.decoding_window(1));
    else
        fprintf('\nWarning: starting timepoint not found, starting from beginning of data...\n');
        lims(1) = 1;
    end;
    if ~isempty(find(round(time,3)==dec_args.decoding_window(end),1))
        lims(2) = find(round(time,3)==dec_args.decoding_window(end));
    else
        fprintf('\nWarning: end timepoint not found, decoding until end of data...\n');
        lims(2) = size(data,2);
    end;
    
else
    lims = [1 size(data,2)];
end;

if size(data,2)<lims(2)
    lims(2) = size(data,2);
end;

data = data(:, lims(1):lims(2), :);

if ~isa(data, 'double')
    data = double(data);
end;

results = struct;

%cross-validation indices - convoluted but can be kept constant if need be
if isempty(dec_args.cv_indices)
    cv = cvpartition(labels, 'kfold', svm_par.kfold);
else
    cv = dec_args.cv_indices;
end;

if isa(cv,'cvpartition')
    cv_train = zeros(length(labels),1); cv_test = cv_train;
    for ii = 1:dec_args.kfold
        cv_train(:,ii) = cv.training(ii);
        cv_test(:,ii) = cv.test(ii);
    end;
else
    if size(cv,2)~=dec_args.kfold
        error('Crossval indices must be supplied in indices x folds matrix or logical array.');
    end;
    cv_train = cv;
    cv_test = abs(cv_train-1);
end;

cv_train = logical(cv_train); cv_test = logical(cv_test);
cv_idx = cv_train; %this will be saved for later - prior to trial averaging

%create pseudo-trials if requested; this is done separately for every test set, ensuring independence
if ~isempty(dec_args.pseudo)
    
    fprintf('\nCreating pseudotrials....\r');
    all_data = cell(1,svm_par.kfold); all_labels = cell(1,svm_par.kfold);
    cv_test_tmp = zeros(1,svm_par.kfold); %we have to redo the crossval indices
    
    for ii = 1:5
        [ps_data, ps_labels] = create_pseudotrials(data(:,:,cv_test(:,ii)), labels(cv_test(:,ii)), dec_args.pseudo(1), dec_args.pseudo(2));
        if ndims(ps_data)>3
            ps_data = reshape(ps_data, size(ps_data,1), size(ps_data,2), size(ps_data,3)*size(ps_data,4));
            ps_labels = reshape(ps_labels, size(ps_labels,1)*size(ps_labels,2),1);
        end;
        all_data{ii} = ps_data; all_labels{ii} = ps_labels;
        if ii==1
            cv_test_tmp(1:length(ps_labels) ,ii) = 1;
            idx = 0;
        else
            idx = idx + length(all_labels{ii-1});
            cv_test_tmp(idx+1:idx+length(ps_labels),ii) = 1;
        end;
    end;
    
    data = cat(3,all_data{:}); labels = cat(1,all_labels{:});
    cv_train = abs(cv_test_tmp-1);
    cv_test = cv_test_tmp;
    cv_train = logical(cv_train); cv_test = logical(cv_test);
    clear all_data all_labels cv_test_tmp;
    
end;

allscore = zeros(length(labels), length(chan_idx), floor(size(data,2)/dec_args.window_length)); accuracy = zeros(3,5, length(chan_idx), floor(size(data,2)/dec_args.window_length));

for ii = 1:svm_par.kfold
    
    fprintf('\rDecoding fold %d out of %d', ii, svm_par.kfold);
    
    if dec_args.mnn
        class_id = unique(labels); class1 = find(labels==class_id(1)); class2 = find(labels==class_id(2)); %this needs to be done separately for each condition
        cv_tmp = find(cv_train(:,ii)==1); class1 = cv_tmp(ismember(cv_tmp,class1)); class2 = cv_tmp(ismember(cv_tmp,class2));
        sigma_time = zeros(2,size(data,2), size(data,1), size(data,1));
        for t = 1:size(data,2)
            sigma_time(1,t,:,:) = cov1para(squeeze(data(:,t,class1))');
            sigma_time(2,t,:,:) = cov1para(squeeze(data(:,t,class2))');
        end;
        sigma_time = squeeze(mean(sigma_time,1)); %average across conditions
        sigma_inv = (squeeze(mean(sigma_time,1)))^-0.5;
        for t = 1:size(data,2)
            data(:,t,cv_train(:,ii)) = (squeeze(data(:,t,cv_train(:,ii)))'*sigma_inv)';
            data(:,t,cv_test(:,ii)) = (squeeze(data(:,t,cv_test(:,ii)))'*sigma_inv)';
            
        end;
    end;
    
    tp = 1:dec_args.window_length:size(data,2)-dec_args.window_length+1;
    fprintf('\nRunning searchlight ');
    
    for c = 1:length(chan_idx)
        
        fprintf('%d out of %d', c, length(chan_idx));
        
        for t = 1:length(tp)
            
            kdata = reshape(data(chan_idx{c},tp(t):tp(t)+dec_args.window_length-1,:), length(chan_idx{c})*dec_args.window_length, size(data,3))'; %select time point or time window
            if svm_par.standardize
                kdata = (kdata - repmat(min(kdata(cv_train(:,ii),:), [], 1), size(kdata, 1), 1)) ./ repmat(max(kdata(cv_train(:,ii),:), [], 1) - min(kdata(cv_train(:,ii),:), [], 1), size(kdata, 1), 1);
            end;
            svm_model = train(labels(cv_train(:,ii)), sparse(kdata(cv_train(:,ii),:)), sprintf('-s %d -c %d -q 1', svm_par.solver, svm_par.boxconstraint)); %dual-problem L2 solver with C=1
            [allscore(cv_test(:,ii),c,t), accuracy(:,ii,c,t), ~] = predict(labels(cv_test(:,ii)), sparse(kdata(cv_test(:,ii),:)), svm_model, '-q 1');
            
        end;
        
        fprintf((repmat('\b',1,numel([num2str(c) num2str(length(chan_idx))])+8)));
        
    end;
    
    fprintf([repmat('\b',1,20) 'Done']);
    
end;

%store a bunch of stuff
results.Accuracy = squeeze(mean(accuracy(1,:,:,:),2));
results.AccuracyMSError = squeeze(mean(accuracy(2,:,:,:),2));
results.AccuracyFold = squeeze(accuracy(1,:,:,:));
for c = 1:length(chan_idx)
    for t = 1:length(tp)
        results.Confusion{c,t} = confusionmat(labels,squeeze(allscore(:,c,t)));
        if numel(results.Confusion)>1
            results.Sensitivity(c,t) = results.Confusion{c,t}(1,1)/(sum(results.Confusion{c,t}(1,:))); %TP/allP = TP/(TP+FN)
            results.Specificity(c,t) = results.Confusion{c,t}(2,2)/(sum(results.Confusion{c,t}(2,:))); %TN/allN = TN/(FP+TN)
            PP = results.Confusion{c,t}(1,1)/(sum(results.Confusion{c,t}(:,1))); %positive predictive value: class1
            NP = results.Confusion{c,t}(2,2)/(sum(results.Confusion{c,t}(:,2))); %negative predictive value: class2
            results.Fscore1(c,t) = (2*PP*results.Sensitivity(c,t))/(PP+results.Sensitivity(c,t));
            results.Fscore2(c,t) = (2*NP*results.Specificity(c,t))/(NP+results.Specificity(c,t));
            results.WeightedFscore(c,t) = ((sum(results.Confusion{c,t}(:,1))/sum(results.Confusion{c,t}(:)))*results.Fscore1(c,t)) + ((sum(results.Confusion{c,t}(:,2))/sum(results.Confusion{c,t}(:)))*results.Fscore2(c,t));
        end;
    end;
end;
results.cv_indices = cv_idx; %this can be reused

end