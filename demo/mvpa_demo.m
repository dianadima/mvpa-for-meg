%% add all scripts to path

addpath('..')
mvpa_setup

%% prepare data for decoding

[data, labels] = prepare_binary_data('data_face.mat','data_scrambled.mat');
load('neighbours.mat')

%% time-resolved decoding using all sensors

sens_result1 = time_resolved_kfold(data, labels, 'mnn', false, 'kfold', 5, 'weights', true, ....                          
                      'window_length', 5, 'decoding_window', [-0.1 1], 'time', neighbours(1).time);

sens_result2 = time_resolved_kfold(data, labels, 'mnn', true, 'pseudo', [3 100], 'kfold', 5, 'weights', true,...
                      'window_length', 5, 'decoding_window', [-0.1 1], 'time', neighbours(1).time);
                  
if ~exist('results','dir'), mkdir ('results'); end
save(fullfile('results','sens_result.mat'), 'sens_result*')

%% plotting over time

%get time axis for decoding window and save indices to timepoints at -0.1 and 0.1 s
time = neighbours(1).time; 
tpstart = find(round(time,3)==-0.1); 
dec_time = neighbours(1).time(tpstart:end); %time axis for decoding window
time = dec_time(1:5:end); %downsampled time axis for decoding window
tp100 = find(round(time,3)==0.1); %100 ms time point

%calculate the standard error of the mean using fold-wise accuracy
sem1 = std(sens_result1.AccuracyFold,1)/sqrt(5);
sem2 = std(sens_result2.AccuracyFold,1)/sqrt(5);

%plot time-resolved accuracy for both analyses in one plot
figure('color','w')
plot_time_results(sens_result1.Accuracy, sem1,'time',time,'smooth',5, 'legend', 'No MNN/avg','ylim',[35 105]); hold on
plot_time_results(sens_result2.Accuracy, sem2,'time',time,'smooth',5, 'legend', 'MNN+avg', 'color', 'b','ylim',[35 105])

%save figure
if ~exist('figures',dir), mkdir('figures'); end
saveas(gcf, fullfile('figures','time-resolved-accuracy.png'))

%% temporal generalization

tg_accuracy1 = temporal_generalization(data,labels,'mnn',false, 'decoding_window', [-0.1 1], 'time', neighbours(1).time); %temporal generalization without MNN/averaging
tg_accuracy2 = temporal_generalization(data,labels,'mnn',true,'pseudo', [3 100], 'decoding_window', [-0.1 1], 'time', neighbours(1).time); %with MNN/averaging
save(fullfile('results','tg_accuracy.mat'),'tg_accuracy*')

%we'll highlight clusters (size>=5) with accuracies over 85 
mask1 = logical(tg_accuracy1>85); %create mask
mask2 = logical(tg_accuracy2>85);

colorlim = [min(min([tg_accuracy1;tg_accuracy2])) max(max([tg_accuracy1;tg_accuracy2]))]; %common color axis

figure('color','w')
subplot(1,2,1)
plot_temporal_generalization(tg_accuracy1,'time',dec_time, 'mask',mask1, 'clustersize',5,'colormap', 'parula','colorbar_label',[],'colorlim', colorlim, 'title', 'No MNN/Ave')
subplot(1,2,2)
plot_temporal_generalization(tg_accuracy2,'time',dec_time, 'mask',mask2, 'clustersize',5,'colormap', 'parula','colorlim', colorlim, 'title', 'MNN + Ave')
saveas(gcf, fullfile('figures','tg-accuracy.png'))

%% plotting patterns based on classifier weights (cf. Haufe et al, 2014)

%plot averaged weight-based patterns from 100 ms onwards
ave_weights1 = mean(sens_result1.WeightPatternsNorm(:,tp100:end),2);
ave_weights2 = mean(sens_result2.WeightPatternsNorm(:,tp100:end),2);

figure;
subplot(1,2,1)
plot_sensor_results(ave_weights1, neighbours, 'window_length', 1, 'time_label', false, 'colorlim', [0 1], 'colormap', 'parula'); colorbar; title('No MNN/Ave', 'FontWeight', 'normal')
subplot(1,2,2)
plot_sensor_results(ave_weights2, neighbours, 'window_length', 1, 'time_label', false, 'colorlim', [0 1], 'colormap', 'parula'); colorbar; title('MNN + Ave', 'FontWeight', 'normal')
saveas(gcf, fullfile('figures','average-weights.png'))

%time-resolved plots (only for MNN/averaged data)
plot_sensor_results(sens_result2.WeightPatternsNorm, neighbours, 'time', time, 'window_length', 0.1, 'colorlim', [0 1], 'colormap', 'parula')
saveas(gcf, fullfile('figures','weights2.png'))
movie_sensor_results(sens_result2.WeightPatternsNorm, neighbours, fullfile('figures','weights2.avi'), 'colorlim', [0 1], 'colormap', 'parula', 'result_type', 'Weight patterns')

%% searchlight decoding

%create separate training & test datasets and use a speedier (but less exhaustive) approach
load('data_face.mat'); load('data_scrambled.mat');
[train_data, train_labels] = prepare_binary_data(data_face(:,:,1:45), data_scrambled(:,:,46:90)); %divide data into halves
[test_data, test_labels] = prepare_binary_data(data_face(:,:,46:90), data_scrambled(:,:,1:45));

%searchlight decoding with longer time windows for speed (without and with MNN/trial averaging)
[srcl_accuracy1, srcl_Fscore1] = searchlight_holdout(train_data, train_labels, test_data, test_labels, neighbours, 'window_length', 10, 'mnn', false);
[srcl_accuracy2, srcl_Fscore2] = searchlight_holdout(train_data, train_labels, test_data, test_labels, neighbours, 'window_length', 10, 'mnn', true, 'pseudo', [3 100]);
save(fullfile('results','srcl_sens_result.mat'),'srcl*')

%the time axis has been downsampled further
srcl_time = neighbours(1).time(1:10:end); %get time axis in steps of 10 samples
tp100 = find(round(srcl_time,3)==0.1); %find the 100 ms time point

%plot movies showing searchlight accuracies over time
colorlim = [min(min([srcl_accuracy1;srcl_accuracy2])) max(max([srcl_accuracy1;srcl_accuracy2]))]; %get common color axis
movie_sensor_results(srcl_accuracy1, neighbours, fullfile('figures','srcl_accuracy1.avi'), 'colormap', 'parula', 'colorlim', colorlim) 
movie_sensor_results(srcl_accuracy2, neighbours, fullfile('figures','srcl_accuracy2.avi'), 'colormap', 'parula', 'colorlim', colorlim)

% plot average searchlight accuracy after 100 ms for both analyses
mean_acc1 = mean(srcl_accuracy1(:,tp100:end),2);
mean_acc2 = mean(srcl_accuracy2(:,tp100:end),2);

colorlim = [min([mean_acc1;mean_acc2]) max([mean_acc1;mean_acc2])];

figure('color','w')
subplot(1,2,1)
plot_sensor_results(mean_acc1, neighbours, 'colormap', 'parula', 'colorlim', colorlim, 'window_length', 1, 'time_label', false); colorbar; title('No MNN/Ave', 'FontWeight', 'normal')
subplot(1,2,2)
plot_sensor_results(mean_acc2, neighbours, 'colormap', 'parula', 'colorlim', colorlim, 'window_length', 1, 'time_label', false); colorbar; title('MNN + Ave', 'FontWeight', 'normal')
saveas(gcf,fullfile('figures','average-srcl-accuracy.png'))

