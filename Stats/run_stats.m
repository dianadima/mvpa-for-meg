function [stats] = run_stats(accuracy, varargin)
% Run sign-shuffling permutation testing (one-tailed) on accuracy measures, correcting for multiple comparisons across time/space/both.
%
% Inputs: accuracy can be: subjects x time, subjects x space x time, subjects x time x time
% NOTE: If data is subjects by space by time, make sure time is last dimension!
%
% Optional inputs: method -- default omnibus, method for correcting for multiple comparisons: cluster, omnibus, fdr, mixed_of (omnibus/fdr along dimensions 2 and 3), mixed_fo (fdr and omnibus along dim 2 and 3 respectively)
%                  spatial_def -- must be provided for spatial cluster correction (sensor-space neighbours structure, or sourcemodel for source-space)
%                  alpha -- default 0.05
%                  clusteralpha -- cluster-defining alpha, default 0.05
%                  num_iterations -- number of sign permutations, default 5000
%                  chance_level -- will be subtracted from accuracy, default 50
%                  statistic -- currently only mean
%                  clusterstatistic -- max, maxsize, maxsum, and wcm, default maxsum
%                  clusterthresh -- common (same threshold across sources/timepoints) or individual (separate thresholds), default common
%
% Output: stats structure containing p-values and other info
%
% DC Dima 2018 (diana.c.dima@gmail.com)

p = inputParser;
addParameter(p, 'method','omnibus'); %cluster, omnibus, fdr, mixed_of, mixed_fo (specify for each dimension)
addParameter(p, 'alpha', 0.05);
addParameter(p, 'clusteralpha', 0.05);
addParameter(p, 'spatial_def', []); %neighbours or sourcemodel for space-resolved data with cluster correction
addParameter(p, 'num_iterations', 5);
addParameter(p, 'chance_level', 50);
addParameter(p, 'statistic', 'mean');
addParameter(p, 'clusterstatistic', 'maxsum');
addParameter(p, 'clusterthresh','individual');
parse(p, varargin{:});
opt = p.Results;

num_it = opt.num_iterations;
alpha = opt.alpha;
spdef = opt.spatial_def;

%get size info
nd = ndims(accuracy); sz = size(accuracy); sz = sz(2:end); if numel(sz)==1, sz = [sz 1]; end; %size without subject dimension
if ~ismember(nd,[2 3])
    error('Accuracy must be a matrix with subjects as first dimension and with maximum 3 dimensions');
end

%create null distribution through sign randomization of accuracies
fprintf('\nEstimating null distribution...')
[r_stat, obs_stat, pval] = randomize_accuracy(accuracy,'num_iterations', num_it, 'chance_level', opt.chance_level, 'statistic', opt.statistic);

switch(opt.method)
    
    case 'omnibus'
        
        fprintf('\nPerforming omnibus correction for multiple comparisons...')
        max_stat = squeeze(max(r_stat,[],2)); if size(max_stat,2)>1, max_stat = squeeze(max(max_stat,[],2)); end
        
        pval = nan(sz); %recalculate pvals using maximum thresholding
        for i = 1:size(accuracy,2)
            for ii = 1:size(accuracy,3)
                pval(i,ii) = (length(find(max_stat>=obs_stat(i,ii)))+1)/(num_it+1);
            end
        end
        stats.pval = pval; stats.method = 'omnibus'; stats.mask = pval<alpha;
        
    case 'fdr'
        
        fprintf('\nPerforming FDR correction for multiple comparisons...')
        [q,mask] = fdr(pval, alpha);
        stats.pval = pval; stats.method = 'fdr'; stats.q = q; stats.mask = mask;
        
    case 'cluster'
        
        fprintf('\nPerforming cluster correction for multiple comparisons...')
        
        %include channel by time, source (grid) by time, and time by time cases
        if ~isempty(spdef)
            
            %we have some spatial info - use FT clusterstat function
            %SPM toolbox needed for spm_bwlabel
            [~,ftpath] = ft_version; private_path = fullfile(ftpath, 'private');
            current_path = pwd;
            cd(private_path)
            
            %common options regardless of data format
            cfg = [];
            cfg.tail = 1;
            cfg.clustertail = 1;
            cfg.feedback = 'yes';
            cfg.clusterstatistic = opt.clusterstatistic;
            cfg.clusteralpha = opt.clusteralpha;
            if strcmp(opt.clusterthresh,'common')
                cfg.clusterthreshold = 'nonparametric_common';
            elseif strcmp(opt.clusterthresh,'individual')
                cfg.clusterthreshold = 'nonparametric_individual';
            end
            cfg.numrandomization = num_it;
            
            if isfield (spdef, 'neighblabel') % channel case
                
                cfgn = [];
                cfgn.neighbours = spdef;
                cfgn.channel = {spdef(:).label};
                conn = channelconnectivity(cfgn);
                
                cfg.dim = size(obs_stat);
                cfg.connectivity = conn;
                
            elseif ~isempty(spdef) && isfield(spdef, 'dim') %sourcemodel case
                
                %create a labelmatrix based on source neighbours
                src = get_source_info;
                connmat = zeros(length(src),length(src));
                for s = 1:length(src)
                    connmat(s,src{s}) = 1;
                end
                
                cfg.dim = size(obs_stat);
                %cfg.origdim = [spdef.dim size(obs_stat,2)];
                %cfg.inside = repmat(spdef.inside,1,size(obs_stat,2));
                cfg.inside = spdef.inside;
                cfg.connectivity = connmat;%conndef(length(spdef.dim),'min'); %can be reshaped to 3D volume
                
            end
            
            rndmap = permute(r_stat,[2 3 1]);
            rndmap = reshape(rndmap,[numel(obs_stat), size(rndmap,3)]);
            clstat = clusterstat(cfg, rndmap, obs_stat(:));
            
            stats.clusterlabelmat = reshape(clstat.posclusterslabelmat, size(obs_stat));
            stats.clusterpvalmat = reshape(clstat.prob,size(obs_stat));
            stats.clusterstat = [clstat.posclusters.clusterstat];
            stats.clusters = cell(1,length(stats.clusterstat));
            for i = 1:length(stats.clusterstat)
                stats.clusters{i} = find(stats.clusterlabelmat==i);
            end
            stats.sigclusters = stats.clusters([clstat.posclusters.prob]<=alpha);
            stats.clusterpvals = [clstat.posclusters.prob];
            stats.randclustermaxdistr = clstat.posdistribution;
            if ismatrix(accuracy)
                stats.sigtime = zeros(1,size(accuracy,2));
                tpidx = cat(1,stats.sigclusters);
                stats.sigtime(tpidx) = 1;
            end
            
            cd(current_path)
            
        else
            
            %temporal/continuous clusters with no spatial structure
            stats = find2Dclusters(obs_stat,r_stat,opt);
            
        end
        
        
    case {'mixed_of', 'mixed_fo'}
        
        fprintf('\nPerforming omnibus and FDR correction for multiple comparisons...')
        if nd==2, error('Mixed methods of MC correction only work when there are at least 2 dimensions, e.g. space x time, or time x time.'); end
        
        %if omnibus correction is desired along the last dimension
        if strcmp(opt.method, 'mixed_fo')
            r_stat = permute(r_stat, [1 3 2]); obs_stat = permute(obs_stat, [1 3 2]);
        end
        
        max_stat = squeeze(max(r_stat,[],2));
        pval = nan(sz); mask = nan(sz); q = nan(1,size(accuracy,2));
        
        for i = 1:size(accuracy,2)
            for ii = 1:size(accuracy,3)
                pval(i,ii) = (length(find(max_stat(:,ii)>=obs_stat(i,ii)))+1)/(num_it+1);
            end
            [q(i), mask(i,:)] = fdr(pval(i,:), alpha); %fdr across 3rd dim
        end
        
        if strcmp(opt.method, 'mixed_fo'), mask = mask'; pval = pval'; end %restore dimension order
        stats.mask = mask; stats.pval = pval;
        
    otherwise
        
        error('Wrong method for MC correction, see help');
        
end

stats.accuracy = accuracy;
stats.method = opt.method;

end




function cl2Dstats = find2Dclusters(obs_stat,r_stat,opt)

num_it = opt.num_iterations;
prc = 100* (1 - opt.clusteralpha); %get cluster-setting percentile
if strcmp(opt.clusterthresh, 'individual')
    thresh = prctile(r_stat,prc,1); %individual threshold setting, as in FT nonparametric_individual option
    o_thresh = squeeze(thresh); if isvector(o_thresh), o_thresh = o_thresh(:); end
    obs_map = double(obs_stat>=o_thresh);
    if ismatrix(r_stat), r_thresh = repmat(thresh,num_it,1); elseif ndims(r_stat)==3, r_thresh = repmat(thresh,num_it,1,1); end
    r_map = double(r_stat>=r_thresh); %ones for values that should go in the clusters
elseif strcmp(opt.clusterthresh,'common') %common threshold across sensors/sources/timepoints
    thresh = prctile(r_stat(:),prc);
    obs_map = double(obs_stat>=thresh);
    r_map = double(r_stat>=thresh);
end
max_r_cls = zeros(1,num_it); %maximal cluster statistic distribution

%here we are simply looking for clusters w/o spatial  structure
conn = conndef(ndims(obs_map),'max'); %use maximal structure as diagonal points belong to same cluster (e.g. temporal generalization)
obs_cls = bwconncomp(obs_map,conn);

for i = 1:num_it
    
    r_map_ = squeeze(r_map(i,:,:));
    r_cls =  bwconncomp(r_map_,conn);
    r_stat_ = squeeze(r_stat(i,:,:));
    if ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'maxsize')
        max_r_cls(i) = max(cellfun(@length,r_cls.PixelIdxList));
    elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'maxsum')
        tmp_cls = nan(1,length(r_cls.PixelIdxList));
        for ii = 1:length(r_cls.PixelIdxList)
            tmp_cls(ii) = sum(r_stat_(r_cls.PixelIdxList{ii}));
        end
        max_r_cls(i) = max(tmp_cls);
    elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'max')
        tmp_cls = nan(1,length(r_cls.PixelIdxList));
        for ii = 1:length(r_cls.PixelIdxList)
            tmp_cls(ii) = max(r_stat_(r_cls.PixelIdxList{ii}));
        end
        max_r_cls(i) = max(tmp_cls);
    elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'wcm')
        tmp_cls = nan(1,length(r_cls.PixelIdxList));
        for ii = 1:length(r_cls.PixelIdxList)
            if strcmp(opt.clusterthresh,'individual')
                t = thresh_(r_cls.PixelIdxList{ii});
            else
                t = thresh;
            end
            tmp_cls(ii) = sum((r_stat_(r_cls.PixelIdxList{ii})-t).^1);
        end
        max_r_cls(i) = max(tmp_cls);
    end
end

%now compare observed with random clusters
if strcmp(opt.clusterstatistic,'maxsize')
    obs_clstat = cellfun(@length, obs_cls.PixelIdxList);
    if ~isempty(obs_clstat)
        cluster_pvals = nan(1,length(obs_clstat));
        for i = 1:length(obs_clstat)
            cluster_pvals(i) = ((sum(max_r_cls>=obs_clstat(i)))+1)/(num_it+1);
        end
    end
elseif strcmp(opt.clusterstatistic,'maxsum')
    cluster_pvals = nan(1,length(obs_cls.PixelIdxList));
    obs_clstat = cluster_pvals;
    for i = 1:length(obs_cls.PixelIdxList)
        obs = obs_stat(obs_cls.PixelIdxList{i});
        obs_clstat(i) = sum(obs);
        cluster_pvals(i) = ((sum(max_r_cls>=obs_clstat(i)))+1)/(num_it+1);
    end
elseif strcmp(opt.clusterstatistic,'max')
    cluster_pvals = nan(1,length(obs_cls.PixelIdxList));
    obs_clstat = cluster_pvals;
    for i = 1:length(obs_cls.PixelIdxList)
        obs = obs_stat(obs_cls.PixelIdxList{i});
        obs_clstat(i) = max(obs);
        cluster_pvals(i) = ((sum(max_r_cls>=obs_clstat(i)))+1)/(num_it+1);
    end
elseif strcmp(opt.clusterstatistic,'wcm')
    cluster_pvals = nan(1,length(obs_cls.PixelIdxList));
    obs_clstat = cluster_pvals;
    for i = 1:length(obs_cls.PixelIdxList)
        obs = obs_stat(obs_cls.PixelIdxList{i});
        if strcmp(opt.clusterthresh,'individual')
            t = thresh_(obs_cls.PixelIdxList{i});
        else
            t = thresh;
        end
        obs_clstat(i) = sum((obs(:)-t(:)).^1);
        cluster_pvals(i) = ((sum(max_r_cls>=obs_clstat(i)))+1)/(num_it+1);
    end
    max_r_cls(i) = max(tmp_cls);
end

%save stuff
fprintf('\nSaving...')
cl2Dstats.clusters = obs_cls.PixelIdxList;
cl2Dstats.clusterstat = obs_clstat;
cl2Dstats.clusterpvals = cluster_pvals;
cl2Dstats.randclusterstatmax = max_r_cls;
cl2Dstats.sigclusters = cl2Dstats.clusters(cl2Dstats.clusterpvals<=opt.alpha);
cl2Dstats.sigpvals = cl2Dstats.clusterpvals(cl2Dstats.clusterpvals<=opt.alpha);

%save an indexing vector for time-resolved data
if isvector(obs_stat)
    cl2Dstats.sigtime = zeros(1,length(obs_stat));
    tpidx = cat(1,cl2Dstats.clusters{cl2Dstats.clusterpvals<=opt.alpha});
    cl2Dstats.sigtime(tpidx) = 1;
end
end


