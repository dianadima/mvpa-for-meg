function [ virtualdata ] = get_lcmv_data( dataset, mri_file, marker, varargin )
% Beamforming of data from one condition, using a LCMV beamformer in Fieldtrip.
% Returns virtual sensor data for all cortical sources.
% Inputs: dataset (path to MEG data); mri_file (path to MRI); marker (trigger name)
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
%        sourcemodel = 'standard_sourcemodel3d10mm'
%
%DC Dima 2017 (diana.c.dima@gmail.com)

addpath(strrep(mfilename('fullpath'), fullfile('beamforming','get_lcmv_data'),''));
opt = preproc.lcmv_args;
list = fieldnames(opt);
p = inputParser;
for i = 1:length(list)
    addParameter(p, list{i}, opt.(list{i}));
end;
addParameter(p,'toilim', [0.5 0.7]); 
addParameter(p,'trlidx', []);    
parse(p, varargin{:});
trialfun = p.Results.trialfun;

[~,ft_path] = ft_version;

%% read in mri

[~,~,mritype] = fileparts(mri_file);
seg_mri = strrep(mri_file, mritype, 'seg.mat');

if strcmp(mritype,'mat')  
    if exist(seg_mri, 'file')
        seg = load(seg_mri);
        mri = load(mri_file);
    elseif ~exist(mri_file, 'file')
        error('Can`t find MRI.')
    else
        mri = load(mri_file);
        %mri = ft_read_mri(raw_mri);
        mri = ft_convert_units(mri, 'mm');
        
        cfg = [];
        seg = ft_volumesegment(cfg, mri); %default - tissue probability map
        seg.transform = mri.transform;
        seg.anatomy = mri.anatomy;
        save(seg_mri, '-struct', 'seg');
    end   
else %assume .mri or .nii
    mri = ft_read_mri(mri_file);
    mri = ft_convert_units(mri, 'mm');
    mri.coordsys = p.Results.coordsys;
    
    cfg = [];
    seg = ft_volumesegment(cfg, mri); %default - tissue probability map
    seg.transform = mri.transform;
    seg.anatomy = mri.anatomy;
    save(seg_mri, '-struct', 'seg');
end

%% read in data
%COND 1
cfg = [];
cfg.dataset = dataset; 
cfg.trialdef.prestim = p.Results.prestim;
cfg.trialdef.poststim = p.Results.poststim;
cfg.trialdef.eventtype = marker;
cfg.trialfun = trialfun;
cfg.trlidx = p.Results.trlidx;
cfg = ft_definetrial(cfg);

%preprocessing options
cfg.channel = 'MEG'; 
cfg.bpfilter = 'yes';
cfg.bpfreq = p.Results.bandpass;
if cfg.bpfreq(1)<0.5 || cfg.bpfreq(2)<10
    cfg.bpfiltord = 3;
end;
data = ft_preprocessing(cfg);

if ~isempty(p.Results.resamplefs)
    cfg = [];
    cfg.resamplefs = p.Results.resamplefs;
    cfg.detrend = 'no';
    data = ft_resampledata(cfg,data);
end;

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
sourcemodel_filename = p.Results.sourcemodel;
sourcemodel_path = fullfile(ft_path,'template','sourcemodel',sourcemodel_filename);
sourcemodel = load(sourcemodel_path);
grid = ft_convert_units(sourcemodel.sourcemodel, 'mm'); 
clear sourcemodel*

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
cfg.trialfun = trialfun; 
cfg.trialdef.eventtype  = marker; 
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
data = ft_preprocessing(cfg);

if ~isempty(p.Results.resamplefs)
    cfg = [];
    cfg.resamplefs = p.Results.resamplefs;
    cfg.detrend = 'no';
    data = ft_resampledata(cfg,data);
end;

% for MNN, calculate error covariance and multiply with the beamformer weights
% Note: work in progress! The beamformer suppresses noise and MNN here 
% might have limited added benefit. 
if p.Results.mnn
    data_ = cat(3, data.trial{:});
    fprintf('\rWhitening data...\r');
    sigma_time = zeros(size(data_,2), size(data_,1), size(data_,1));
    for t = 1:size(data_,2)
        sigma_time(t,:,:) = cov1para(squeeze(data_(:,t,:))');
    end;
    sigma_inv = (squeeze(mean(sigma_time,1)))^-0.5;
    clear data_;
end;

virtualdata = data;
%Create the virtual data
for j = 1 : length(data.trial)
    if p.Results.fixedori
        if p.Results.mnn
            virtualdata.trial{j} = (filters*sigma_inv) * data.trial{j};
        else
            virtualdata.trial{j} = filters * data.trial{j};
        end;
    else
        virtualdata.trial{j} = zeros(size(filters,1),size(data.trial{j},2),size(filters,3));
        for ori = 1:3
            if p.Results.mnn
                virtualdata.trial{j}(:,:,ori) = (filters(:,:,ori)*sigma_inv) * data.trial{j};
            else
                virtualdata.trial{j}(:,:,ori) = filters(:,:,ori) * data.trial{j};
            end;
        end
    end;
end;

for j = 1 : size(virtualdata.trial{1},1)    %to keep FT happy
       virtualdata.label{j,1} = num2str(j);
end

virtualdata = rmfield(virtualdata, 'cfg');


end
