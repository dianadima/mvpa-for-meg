function [ results ] = time_resolved_kfold( data, labels, varargin )
% Performs time-resolved SVM decoding of MEG data, using stratified k-fold cross-validation on each time window, and LibLinear SVM implementation.
% Inputs: data, labels.
% Optional: 'sensor_idx', structure obtained using get_sensor_info - for channel selection. You can also just provide numerical indices, in which case you don't need the structure.
%                        If you want to subselect features on source space data, you need to manually provide numerical indice (i.e. 'channels', [1:1000]).
%          'channels', channel set set (string or cell array of strings; default: 'MEG').
%          'decoding_window' (limits; default: [] - all timepoints). In  sampled time points (OR in seconds - only if you also provide time axis).
%          'window_length' (in sampled time points; default: 1).
%          'time', time axis, if you want to give the decoding window in seconds, you also need to provide a time axis, matching the second dimension of the data).
%             
%           * other possible name-value pairs: SVM settings, svm evaluation metrics (see Documentation).
% Outputs: structure containing classification performance metrics (for each timepoint).


%parse inputs
dec_args = decoding_args;
svm_par = svm_args;
list = [fieldnames(dec_args); fieldnames(svm_par)];
p = inputParser;

for i = 1:length(properties(decoding_args))
    addParameter(p, list{i}, dec_args.(list{i}));
end;
for ii = i+1:length(properties(dec_args))+length(properties(svm_args))
    addParameter(p, list{ii}, svm_par.(list{ii}));
end;
addParameter(p, 'sensor_idx', []);


parse(p, varargin{:});
dec_args = p.Results;
svm_par = rmfield(struct(dec_args), {'window_length','channels','decoding_window', 'time', 'sensor_idx'}); %converted struct will be fed into decoding function
clear p;

%get channel indices and time axis. Numerical channel indices take priority
if ~iscell(dec_args.channels) && ~ischar(dec_args.channels)
    chan_idx = dec_args.channels;
end;
if ~isempty(dec_args.sensor_idx) %a structure was given
    sensor_idx = dec_args.sensor_idx;
    if ~exist('chan_idx', 'var') %there are no numerical indices
        chan_idx = 1:size(data,1); %initialize with entire array
        if ~strcmp (dec_args.channels, 'MEG') %if we need to subselect sensors
            chan = [];
            for i = 1:length(dec_args.channels)
                idx = cellfun('isempty',strfind({sensor_idx.label},dec_args.channels{i}));
                chan = [chan chan_idx(~idx)]; %#ok<AGROW>
            end;
            chan_idx = chan;
        end;
    end;
else
    if ~exist('chan_idx', 'var')
        chan_idx = 1:size(data,1);
    end;
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
    fprintf('Warning: time axis does not match dataset size. Replacing with default time axis...');
end;

%subselect data corresponding to decoding window requested
if ~isempty(dec_args.decoding_window)
    if ~isempty(find(round(time,3)==dec_args.decoding_window(1),1))
        lims(1) = find(round(time,3)==dec_args.decoding_window(1));
    else
        fprintf('Warning: starting timepoint not found, starting from beginning of data...\n');
        lims(1) = 1;
    end;
    if ~isempty(find(round(time,3)==dec_args.decoding_window(2),1))
        lims(2) = find(round(time,3)==dec_args.decoding_window(2));
    else
        fprintf('Warning: end timepoint not found, decoding until end of data...\n');
        lims(2) = size(data,2);
    end;
    
else
    lims = [1 size(data,2)];
end;

if size(data,2)<lims(2)
    lims(2) = size(data,2);
end;

data = data(chan_idx, lims(1):lims(2), :);

if ~isa(data, 'double')
    data = double(data);
end;

results = struct;

for icv = 1: svm_par.iterate_cv
    
    %cross-validation indices - convoluted but can be kept constant if need be
    if isempty(dec_args.cv_indices)
        cv = cvpartition(labels, 'kfold', svm_par.kfold);
    else
        if iscell(dec_args.cv_indices)
            cv = dec_args.cv_indices{icv};
        else
            cv = dec_args.cv_indices;
        end;
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
        all_data = cell(1,5); all_labels = cell(1,5);
        cv_test_tmp = zeros(1,5); %we have to redo the crossval indices
       
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
    
    allscore = zeros(length(labels),size(data,2)); accuracy = zeros(3,5, floor(size(data,2)/dec_args.window_length));
    fprintf(['\rDecoding fold ' repmat(' ',1,numel([num2str(ii) num2str(length(svm_par.kfold))])+8)]);
    
    for ii = 1:svm_par.kfold

        fprintf((repmat('\b',1,numel([num2str(ii) num2str(length(svm_par.kfold))])+8)));
                
        if dec_args.mnn
            sigma_time = zeros(size(data,2), size(data,1), size(data,1));
            for t = 1:size(data,2)
                sigma_time(t,:,:) = cov1para(squeeze(data(:,t,cv_train(:,ii)))');
            end;
            sigma_inv = (squeeze(mean(sigma_time,1)))^-0.5;
            for t = 1:size(data,2)
                data(:,t,cv_train(:,ii)) = (squeeze(data(:,t,cv_train(:,ii)))'*sigma_inv)';
                data(:,t,cv_test(:,ii)) = (squeeze(data(:,t,cv_test(:,ii)))'*sigma_inv)';
                
            end;
        end;
        
        tp = 1:dec_args.window_length:size(data,2)-dec_args.window_length+1;
        fprintf('%d out of %d', ii, svm_par.kfold); 
        
        for t = 1:length(tp)
            
            kdata = reshape(data(:,tp(t):tp(t)+dec_args.window_length-1,:), size(data,1)*dec_args.window_length, size(data,3))'; %select time point or time window
            if svm_par.standardize
                kdata_ = (kdata - repmat(min(kdata(cv_train(:,ii),:), [], 1), size(kdata, 1), 1)) ./ repmat(max(kdata(cv_train(:,ii),:), [], 1) - min(kdata(cv_train(:,ii),:), [], 1), size(kdata, 1), 1);
            end;
            svm_model = train(labels(cv_train(:,ii)), sparse(kdata_(cv_train(:,ii),:)), sprintf('-s %d -c %d -q 1', svm_par.solver, svm_par.boxconstraint)); %dual-problem L2 solver with C=1
            [allscore(cv_test(:,ii),t), accuracy(:,ii,t), ~] = predict(labels(cv_test(:,ii)), sparse(kdata_(cv_test(:,ii),:)), svm_model, '-q 1');
            
        end;
        
    end;
    
    %store a bunch of stuff
    results.Accuracy(icv,:) = mean(accuracy(1,:,:),2);
    results.AccuracyMSError(icv,:) = mean(accuracy(2,:,:),2);
    results.AccuracyFold(icv,:,:) = accuracy(1,:,:);
    for t = 1:length(tp)
        results.Confusion{icv,t} = confusionmat(labels,allscore(:,t));
        if numel(results.Confusion{icv,t})>1
            results.Sensitivity(icv,t) = results.Confusion{icv,t}(1,1)/(sum(results.Confusion{icv,t}(1,:))); %TP/allP = TP/(TP+FN)
            results.Specificity(icv,t) = results.Confusion{icv,t}(2,2)/(sum(results.Confusion{icv,t}(2,:))); %TN/allN = TN/(FP+TN)
            PP = results.Confusion{icv,t}(1,1)/(sum(results.Confusion{icv,t}(:,1))); %positive predictive value: class1
            NP = results.Confusion{icv,t}(2,2)/(sum(results.Confusion{icv,t}(:,2))); %negative predictive value: class2
            results.Fscore1(icv,t) = (2*PP*results.Sensitivity(icv,t))/(PP+results.Sensitivity(icv,t));
            results.Fscore2(icv,t) = (2*NP*results.Specificity(icv,t))/(NP+results.Specificity(icv,t));
            results.WeightedFscore(icv,t) = ((sum(results.Confusion{icv,t}(:,1))/sum(results.Confusion{icv,t}(:)))*results.Fscore1(icv,t)) + ((sum(results.Confusion{icv,t}(:,2))/sum(results.Confusion{icv,t}(:)))*results.Fscore2(icv,t));
        end;
    end;
    results.cv_indices(icv,:,:) = cv_idx; %this can be reused
    
    
    %calculate weights and compute activation patterns as per Haufe (2014)
    if svm_par.weights
        if dec_args.mnn
            sigma_time = zeros(size(data,2), size(data,1), size(data,1));
            for t = 1:size(data,2)
                sigma_time(t,:,:) = cov1para(squeeze(data(:,t,:))');
            end;
            sigma_inv = (squeeze(mean(sigma_time,1)))^-0.5;
            for t = 1:size(data,2)
                data(:,t,:) = (squeeze(data(:,t,:))'*sigma_inv)';
            end;
        end;
        for t = 1:length(tp)
            kdata = reshape(data(:,tp(t):tp(t)+dec_args.window_length-1,:), size(data,1)*dec_args.window_length, size(data,3))'; %select time point or time window
            if svm_par.standardize
                kdata_ = (kdata - repmat(min(kdata, [], 1), size(kdata, 1), 1)) ./ repmat(max(kdata, [], 1) - min(kdata, [], 1), size(kdata, 1), 1);
            end;
            svm_model = train(labels, sparse(kdata_), sprintf('-s %d -c %d -q 1', svm_par.solver, svm_par.boxconstraint));
            results.Weights(:,t) = svm_model.w;
            results.WeightPatterns(:,t) = abs(cov(kdata_)*results.Weights(:,t)/cov(kdata_*results.Weights(:,t)));
        end;
        results.WeightPatternsNorm = (results.WeightPatterns-min(results.WeightPatterns(:)))/(max(results.WeightPatterns(:))-min(results.WeightPatterns(:)));
    end;
    
    
end
