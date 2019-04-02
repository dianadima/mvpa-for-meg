function [ virtualdata1, virtualdata2 ] = get_lcmv_data_aal( dataset, mri_file, marker1, marker2, varargin )
% Beamforming of data from two conditions with a common filter, using a LCMV beamformer in Fieldtrip.
% Returns virtual sensor data for 90 peak voxels in AAL regions within required frequency band.
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
%        sourcemodel = 'standard_sourcemodel3d10mm'
%
%DC Dima 2017 (diana.c.dima@gmail.com)

addpath(strrep(mfilename('fullpath'), fullfile('beamforming','get_lcmv_data_aal'),''));
opt = preproc.lcmv_args;
list = fieldnames(opt);
p = inputParser;
for i = 1:length(list)
    addParameter(p, list{i}, opt.(list{i}));
end;
addParameter(p,'toilim', [0.5 0.7]); 
addParameter(p,'trlidx', []);    
parse(p, varargin{:});

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
cfg.dataset = dataset; 
cfg.trialdef.prestim = p.Results.prestim;
cfg.trialdef.poststim = p.Results.poststim;
cfg.trialdef.eventtype = marker1;
cfg.trialfun = p.Results.trialfun;
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
cfg.trialfun = p.Results.trialfun;
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
cfg.lcmv.fixedori = 'yes';
cfg.lcmv.projectnoise='yes';
cfg.lcmv.lambda = '5%';
source = ft_sourceanalysis(cfg, data_tl);
source.pos = sourcemodel.pos; %replace coordinates with template ones

%normalize filters if requested (Hillebrand et al., 2012)
filters = cell2mat(source.avg.filter(source.inside == 1));
if p.Results.normalize
    for i = 1 : size(filters,1) %Uses the vector norm for weights normalisation (S.M.)
       filters(i,:) = filters(i,:) ./ norm(filters(i,:));
    end
end

%now the data to multiply with the filters
cfg = []; 
cfg.demean = 'yes'; 
cfg.baselinewindow = p.Results.baseline; %baselining
cfg.dataset = dataset;
cfg.trialfun = p.Results.trialfun; 
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
data1 = rmfield(data1, 'cfg');

cfg = []; 
cfg.demean = 'yes'; 
cfg.baselinewindow = p.Results.baseline; %baselining
cfg.dataset = dataset;
cfg.trialfun = p.Results.trialfun; 
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
data2 = rmfield(data2, 'cfg');

if ~isempty(p.Results.resamplefs)
    cfg = [];
    cfg.resamplefs = p.Results.resamplefs;
    cfg.detrend = 'no';
    data1 = ft_resampledata(cfg,data1);
    data2 = ft_resampledata(cfg,data2);
end;

%Now get the AAL peaks for combined data
[aal_idx, roi_labels] = get_source_info('source_selection', 'aal90');
aal_max = zeros(1,length(aal_idx));
 for roi = 1:length(aal_idx)
     fprintf('\nFinding peak voxel in ROI%d out of 90...\n', roi);
     vdata = data;
     for i = 1 : length(vdata.trial)
         vdata.trial{i} = filters(aal_idx{roi},:) * vdata.trial{i};
     end;
     vdata.label = cell(1,size(vdata.trial{1},1));
     for i = 1:length(vdata.label)
         vdata.label{i} = num2str(i);
     end;
     
     cfg = [];
     cfg.foilim = p.Results.bandpass;
     cfg.method = 'mtmfft';
     cfg.taper = 'hanning';
     freq = ft_freqanalysis(cfg,vdata);
     
     [~,index] = max(mean(freq.powspctrm,2));
     aal_max(roi) = aal_idx{roi}(index);
 end;

 virtualdata1 = data1;
%Create the virtual data
for j = 1 : length(data1.trial)    
       virtualdata1.trial{j} = filters(aal_max,:) * data1.trial{j};
end
virtualdata1.label = roi_labels;

virtualdata2 = data2;
%Create the virtual data
for j = 1 : length(data2.trial)    
       virtualdata2.trial{j} = filters(aal_max,:) * data2.trial{j};
end
virtualdata2.label = roi_labels;

end
