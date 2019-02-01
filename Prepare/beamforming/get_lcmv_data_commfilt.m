function [ virtualdata1, virtualdata2 ] = get_lcmv_data_commfilt( dataset, mri_file, marker1, marker2, varargin )
% Beamforming of data from two conditions with a common filter, using a LCMV beamformer in Fieldtrip.
% Returns virtual sensor data for all cortical sources.
% Inputs: dataset (path to MEG data); mri_file (path to MRI); marker1 and marker2 (trigger names)
% Optional name-value inputs with their defaults:
%        normalize = true;
%        marker = 'onset';
%        trialfun = 'ft_trialfun_general';
%        bandpass = [0.1 100];
%        prestim = 0.5;
%        poststim = 0.7;
%        baseline = [-0.5 0];
%        resamplefs = 300;
%        fixedori = true;
%        plot = 0;
%        mnn = false; (multivariate noise normalization - Work In Progress)
%
%DC Dima 2017 (diana.c.dima@gmail.com)

addpath(strrep(mfilename('fullpath'), fullfile('beamforming','get_lcmv_data_commfilt'),''));

opt = preproc.lcmv_args;
list = fieldnames(opt);
p = inputParser;
for i = 1:length(list)
    addParameter(p, list{i}, opt.(list{i}));
end;
addParameter(p,'toilim', [0.5 0.7]); 
addParameter(p,'trlidx', []);    
parse(p, varargin{:});

%can read the datasets in using different trial functions (random but can be useful)
if iscell(p.Results.trialfun)
    trialfun1 = p.Results.trialfun{1}; trialfun2 = p.Results.trialfun{2};
else
    trialfun1 = p.Results.trialfun; trialfun2 = trialfun1;
end;

%% read in mri

%[path,~] = fileparts(dataset); 
raw_mri = strrep(mri_file, '.mat', '_realigned.mat');

if exist(mri_file, 'file')
    seg = load(mri_file);
    mri = load(raw_mri);
else
    mri = load(raw_mri);
    %mri = ft_read_mri(raw_mri);
    mri = ft_convert_units(mri, 'mm');
    
    cfg = [];
    seg = ft_volumesegment(cfg, mri); %default - tissue probability map
    seg.transform = mri.transform;
    seg.anatomy = mri.anatomy; 
    save(seg_mri, '-struct', 'seg');
end;

%% read in data
%COND 1
cfg = [];
cfg.dataset = dataset; 
cfg.trialdef.prestim = p.Results.prestim;
cfg.trialdef.poststim = p.Results.poststim;
cfg.trialdef.eventtype = marker1;
cfg.trialfun = trialfun1;
cfg.trlidx = p.Results.trlidx;
cfg = ft_definetrial(cfg);

%preprocessing options
cfg.channel = 'MEG'; 
cfg.bpfilter = 'yes';
cfg.bpfreq = p.Results.bandpass;
if cfg.bpfreq(1)<0.5 || cfg.bpfreq(2)<10
    cfg.bpfiltord = 3;
end;
data1 = ft_preprocessing(cfg);

%COND 2
cfg = [];
cfg.dataset = dataset; 
cfg.trialdef.prestim = p.Results.prestim;
cfg.trialdef.poststim = p.Results.poststim;
cfg.trialdef.eventtype = marker2;
cfg.trialfun = trialfun2;
cfg.trlidx = p.Results.trlidx;
cfg = ft_definetrial(cfg);

%preprocessing options
cfg.channel = 'MEG'; 
cfg.bpfilter = 'yes';
cfg.bpfreq = p.Results.bandpass;
if cfg.bpfreq(1)<0.5 || cfg.bpfreq(2)<10
    cfg.bpfiltord = 3;
end;
data2 = ft_preprocessing(cfg);

if ~isempty(p.Results.resamplefs)
    cfg = [];
    cfg.resamplefs = p.Results.resamplefs;
    cfg.detrend = 'no';
    data1 = ft_resampledata(cfg,data1);
    data2 = ft_resampledata(cfg,data2);
end;

%concatenate for common filter
data = ft_appenddata([], data1, data2); 

%timelock analysis
cfg = [];
cfg.channel = 'MEG';
cfg.covariance ='yes';
cfg.removemean = 'no';
data_tl = ft_timelockanalysis(cfg,data);
clear data;    

% reads gradiometer info
hdr = ft_read_header(dataset);

% make head model in indiv space
cfg = [];
cfg.method = 'singleshell'; %localspheres does not seem to align well
cfg.feedback = 'no';
cfg.grad = hdr.grad;
hdm = ft_prepare_headmodel(cfg, seg); %headmodel based on segmented mri
hdm = ft_convert_units(hdm, 'mm');%this is in mm

% get template source model so that individual grids align to MNI
load('/cubric/software/MEG/fieldtrip-20161011/template/sourcemodel/standard_sourcemodel3d10mm');
grid = ft_convert_units(sourcemodel, 'mm'); %#ok<NODEF>

%create sourcemodel by warping individual grid to template grid
cfg = [];
cfg.grad = hdr.grad;
cfg.grid.warpmni= 'yes';
cfg.grid.template = grid;
cfg.grid.nonlinear = 'yes';
cfg.grid.unit = 'mm';
cfg.mri = mri; %sourcemodel based on base mri
sourcemodel = ft_prepare_sourcemodel(cfg); %also mm
clear mri segmentedmri;

% check headmodel
if p.Results.plot==1
     hdm_cm = ft_convert_units(hdm, 'cm'); sourcemodel_cm = ft_convert_units(sourcemodel,'cm');
     figure; hold on     % plot all objects in one figure
     ft_plot_vol(hdm_cm, 'facecolor','skin','edgecolor', 'none','facealpha',0.4)
     ft_plot_mesh(sourcemodel_cm.pos(sourcemodel_cm.inside,:));
     ft_plot_sens(grad);
     hold off
end

% compute leadfield
cfg = [];
cfg.headmodel = hdm;
cfg.reducerank = 2;
cfg.grid = sourcemodel;
cfg.grid.inside = sourcemodel.inside;
cfg.grad = hdr.grad;
cfg.channel = {'MEG'};
lf = ft_prepare_leadfield(cfg); %it's fine to have sensors in cm as this converts between units

%create global spatial filter
cfg = [];
cfg.method = 'lcmv';
cfg.grid = lf;
cfg.channel = {'MEG'};
cfg.headmodel = hdm;
cfg.grad = hdr.grad;
cfg.lcmv.keepfilter = 'yes';
if p.Results.fixedori
    cfg.lcmv.fixedori = 'yes';
else
    cfg.lcmv.fixedori = 'no';
end;
cfg.lcmv.projectnoise='yes';
cfg.lcmv.lambda = '5%';
source = ft_sourceanalysis(cfg, data_tl);
source.pos = sourcemodel.pos; %replace coordinates with template ones

%normalize filters if requested (Hillebrand et al., 2012)
if p.Results.fixedori
    filters = cell2mat(source.avg.filter(source.inside == 1));
    if p.Results.normalize
        for i = 1 : size(filters,1) %Uses the vector norm for weights normalisation (S.M.)
            filters(i,:) = filters(i,:) ./ norm(filters(i,:));
        end
    end
else
    filters = source.avg.filter(source.inside == 1);
    filters = cat(3,filters{:});
    filters = permute(filters, [3 2 1]); %source x time x orientation
    if p.Results.normalize
        for ori = 1:size(filters,3)
            for i = 1 : size(filters,1) %Uses the vector norm for weights normalisation (S.M.)
                filters(i,:,ori) = filters(i,:,ori) ./ norm(filters(i,:,ori));
            end
        end;
    end
end;


%now the data to multiply with the filters
cfg = []; 
cfg.demean = 'yes'; 
cfg.baselinewindow = p.Results.baseline; %baselining
cfg.dataset = dataset;
cfg.trialfun = trialfun1; 
cfg.trialdef.eventtype  = marker1; 
cfg.trialdef.prestim = p.Results.toilim(1);
cfg.trialdef.poststim = p.Results.toilim(2);
cfg.trlidx = p.Results.trlidx;
cfg.channel = 'MEG';
cfg.bpfilter = 'yes';
cfg.bpfreq = p.Results.bandpass;
if cfg.bpfreq(1)<0.5 || cfg.bpfreq(2)<10
    cfg.bpfiltord = 3;
end;
cfg = ft_definetrial(cfg);  
data1 = ft_preprocessing(cfg);


cfg = []; 
cfg.demean = 'yes'; 
cfg.baselinewindow = p.Results.baseline; %baselining
cfg.dataset = dataset;
cfg.trialfun = trialfun2; 
cfg.trialdef.eventtype  = marker2; 
cfg.trialdef.prestim = p.Results.toilim(1);
cfg.trialdef.poststim = p.Results.toilim(2);
cfg.trlidx = p.Results.trlidx;
cfg.channel = 'MEG';
cfg.bpfilter = 'yes';
cfg.bpfreq = p.Results.bandpass;
if cfg.bpfreq(1)<0.5 || cfg.bpfreq(2)<10
    cfg.bpfiltord = 3;
end;
cfg = ft_definetrial(cfg);  
data2 = ft_preprocessing(cfg);


if ~isempty(p.Results.resamplefs)
    cfg = [];
    cfg.resamplefs = p.Results.resamplefs;
    cfg.detrend = 'no';
    data1 = ft_resampledata(cfg,data1);
    data2 = ft_resampledata(cfg,data2);
end;

% for MNN, calculate error covariance based on data from both conditions
% Note: work in progress! This procedure might remove variance of
% interest between two conditions, but it is akin to creating a common filter 
% for beamforming. Furthermore, the beamformer suppresses noise and MNN here 
% might have limited added benefit. 
if p.Results.mnn
    data = ft_appenddata([], data1, data2); 
    data = cat(3, data.trial{:});
    fprintf('\rWhitening data...\r');
    sigma_time = zeros(size(data,2), size(data,1), size(data,1));
    for t = 1:size(data,2)
        sigma_time(t,:,:) = cov1para(squeeze(data(:,t,:))');
    end;
    sigma_inv = (squeeze(mean(sigma_time,1)))^-0.5;
    clear data;
end;

virtualdata1 = data1;
%Create the virtual data
for j = 1 : length(data1.trial)
    if p.Results.fixedori
        if p.Results.mnn
            virtualdata1.trial{j} = (filters*sigma_inv) * data1.trial{j};
        else
            virtualdata1.trial{j} = filters * data1.trial{j};
        end;
    else
        virtualdata1.trial{j} = zeros(size(filters,1),size(data1.trial{j},2),size(filters,3));
        for ori = 1:3
            if p.Results.mnn
                virtualdata1.trial{j}(:,:,ori) = (filters(:,:,ori)*sigma_inv) * data1.trial{j};
            else
                virtualdata1.trial{j}(:,:,ori) = filters(:,:,ori) * data1.trial{j};
            end;
        end
    end;
end;

for j = 1 : size(virtualdata1.trial{1},1)    %to keep FT happy
       virtualdata1.label{j,1} = num2str(j);
end

virtualdata2 = data2;
%Create the virtual data
for j = 1 : length(data2.trial)
    if p.Results.fixedori
        if p.Results.mnn
            virtualdata2.trial{j} = (filters*sigma_inv) * data2.trial{j};
        else
            virtualdata2.trial{j} = filters * data2.trial{j};
        end;
    else
        virtualdata2.trial{j} = zeros(size(filters,1),size(data2.trial{j},2),size(filters,3));
        for ori = 1:3
            if p.Results.mnn
                virtualdata2.trial{j}(:,:,ori) = (filters(:,:,ori)*sigma_inv) * data2.trial{j};
            else
                virtualdata2.trial{j}(:,:,ori) = filters(:,:,ori) * data2.trial{j};
            end;
        end
    end;
end;

for j = 1 : size(virtualdata2.trial{1},1)    %to keep FT happy
       virtualdata2.label{j,1} = num2str(j);
end

virtualdata1 = rmfield(virtualdata1, 'cfg');
virtualdata2 = rmfield(virtualdata2, 'cfg');

end

