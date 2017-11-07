function [ virtualdata ] = get_lcmv_data( dataset, mri_file, varargin )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

opt = lcmv_args;
list = fieldnames(opt);
parsedargs = inputParser;
for i = 1:length(list)
    addParameter(parsedargs, list{i}, opt.(list{i}));
end;
    
parse(parsedargs, varargin{:});

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
cfg.dataset = dataset; 
cfg.trialdef.prestim = parsedargs.Results.prestim;
cfg.trialdef.poststim = parsedargs.Results.poststim;
cfg.trialdef.eventtype = parsedargs.Results.marker;
cfg.trialfun = 'ft_trialfun_NoBadTrials';
cfg = ft_definetrial(cfg);

%preprocessing options
cfg.channel = 'MEG'; 
cfg.bpfilter = 'yes';
cfg.bpfreq = parsedargs.Results.bandpass;
if cfg.bpfreq(1)<0.5 || cfg.bpfreq(2)<10
    cfg.bpfiltord = 3;
end;
data = ft_preprocessing(cfg);

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
if parsedargs.Results.plot==1
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

%now the data to multiply with the filters
cfg = []; 
cfg.demean = 'yes'; 
cfg.baselinewindow = parsedargs.Results.baseline; %baselining
cfg.dataset = dataset;
cfg.trialfun = 'ft_trialfun_NoBadTrials'; 
cfg.trialdef.eventtype  = parsedargs.Results.marker;; 
cfg.trialdef.prestim = parsedargs.Results.prestim;
cfg.trialdef.poststim = parsedargs.Results.poststim;
cfg.channel = 'MEG';
cfg.bpfilter = 'yes';
cfg.bpfreq = parsedargs.Results.bandpass;
if cfg.bpfreq(1)<0.5 || cfg.bpfreq(2)<10
    cfg.bpfiltord = 3;
end;
cfg.precision = 'single';
cfg = ft_definetrial(cfg);  
data = ft_preprocessing(cfg);
data = rmfield(data, 'cfg');

%normalize filters if requested (Hillebrand et al., 2012)
filters = cell2mat(source.avg.filter(source.inside == 1));
if parsedargs.Results.normalize
    for i = 1 : size(filters,1) %Uses the vector norm for weights normalisation (S.Muthukumaraswamy)
       filters(i,:) = filters(i,:) ./ norm(filters(i,:));
    end
end

virtualdata = data;
%Create the virtual data
for j = 1 : length(data.trial)    
       virtualdata.trial{j} = filters * data.trial{j};
end

for j = 1 : size(virtualdata.trial{1},1)    %to keep FT happy
       virtualdata.label{j,1} = num2str(j);
end

%finally downsample if requested
if ~isempty(parsedargs.Results.resamplefs)
    cfg = [];
    cfg.resamplefs = parsedargs.Results.resamplefs;
    cfg.detrend = 'yes';
    cfg.resamplemethod = 'downsample';
    virtualdata = ft_resampledata(cfg,virtualdata);
end;



end

