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
%                  clusterstatistic -- max, maxsize, maxsum, and wcm, default maxsize
%
% Output: stats structure containing p-values and other info
%
% DC Dima 2018 (diana.c.dima@gmail.com)

p = inputParser;
addParameter(p, 'method','omnibus'); %cluster, omnibus, fdr, mixed_of, mixed_fo (specify for each dimension)
addParameter(p, 'alpha', 0.05);
addParameter(p, 'clusteralpha', 0.05);
addParameter(p, 'spatial_def', []); %neighbours or sourcemodel for space-resolved data with cluster correction
addParameter(p, 'num_iterations', 5000);
addParameter(p, 'chance_level', 50);
addParameter(p, 'statistic', 'mean');
addParameter(p, 'clusterstatistic', 'maxsum'); % or 'maxsize'
parse(p, varargin{:});
opt = p.Results;

num_it = opt.num_iterations;
alpha = opt.alpha;

%get size info
nd = ndims(accuracy); sz = size(accuracy); sz = sz(2:end); if numel(sz)==1, sz = [sz 1]; end; %size without subject dimension
if ~ismember(nd,[2 3])
    error('Accuracy must be a matrix with subjects as first dimension and with maximum 3 dimensions');
end

%create null distribution through sign randomization of accuracies
[r_stat, obs_stat, pval] = randomize_accuracy(accuracy,'num_iterations', num_it, 'chance_level', opt.chance_level, 'statistic', opt.statistic);

switch(opt.method)
    
    case 'omnibus'
        
        max_stat = squeeze(max(r_stat,[],2)); if size(max_stat,2)>1, max_stat = squeeze(max(max_stat,[],2)); end
        
        pval = nan(sz); %recalculate pvals using maximum thresholding
        for i = 1:size(accuracy,2)
            for ii = 1:size(accuracy,3)
                pval(i,ii) = (length(find(max_stat>=obs_stat(i,ii)))+1)/(num_it+1);
            end
        end
        stats.pval = pval; stats.method = 'omnibus'; stats.mask = pval<alpha;
        
    case 'fdr'
        
        [q,mask] = fdr(pval, alpha);
        stats.pval = pval; stats.method = 'fdr'; stats.q = q; stats.mask = mask;
        
    case 'cluster'
        
        %look for positive clusters in observed and random data: one-tailed
        prc = 100* (1 - opt.clusteralpha); %get cluster-setting percentile
        thresh = prctile(r_stat(:),prc);
        obs_map = double(obs_stat>=thresh); 
        r_map = double(r_stat>=thresh); %ones for values that should go in the clusters
        max_r_cls = zeros(1,num_it); %maximal cluster statistic distribution
        
        %include channel, source (grid) data, and time by time cases
        if ~isempty(opt.spatial_def) && isfield (opt.spatial_def, 'neighblabel') % channel case
            
            %here we will use some Fieldtrip functions to not reinvent the wheel
            %SPM toolbox needed for spm_bwlabel
            [~,ftpath] = ft_version; private_path = fullfile(ftpath, 'private');
            current_path = pwd;
            
            cfg = [];
            cfg.neighbours = opt.spatial_def;
            cfg.channel = {opt.spatial_def(:).label};
            cd(private_path);
            conn = channelconnectivity(cfg);
            
            cfg = [];
            cfg.dim = size(obs_stat);
            cfg.connectivity = conn;
            cfg.tail = 1; 
            cfg.clustertail = 1;
            cfg.feedback = 'yes';
            cfg.clusterstatistic = opt.clusterstatistic;
            cfg.clusteralpha = opt.clusteralpha;
            cfg.clustercritval = prctile(obs_stat(:),prc);
            cfg.clusterthreshold = 'parametric';
            cfg.numrandomization = size(r_stat,1);
            rndmap = permute(r_stat,[2 3 1]); rndmap = reshape(rndmap,[numel(obs_stat), size(rndmap,3)]);
            clstat = clusterstat(cfg, rndmap, obs_stat(:));
            
            stats.clusterlabels = reshape(clstat.posclusterslabelmat, size(obs_stat));
            stats.clusterpvals = reshape(clstat.prob,size(obs_stat));
            stats.clustersizes = [clstat.posclusters.clusterstat];
            stats.clusters = cell(1,length(stats.clustersizes));
            for i = 1:length(stats.clustersizes), stats.clusters{i} = find(stats.clusterlabels==i); end
            stats.randclustermaxdistr = clstat.posdistribution;
            
            cd(current_path);
            
        else
            
            if ~isempty(opt.spatial_def) && isfield(opt.spatial_def, 'dim') %sourcemodel case
                
                if ~isvector(obs_map), ndim = 1; else ndim = 0; end %space x time case
                
                conn = conndef(length(opt.spatial_def.dim)+ndim, 'max');
                obs_tmp = zeros([opt.spatial_def.dim size(obs_map,2)]);
                for j = 1:size(obs_map,2)
                    obs_tmp(opt.spatial_def.inside,j) = obs_map(:,j);
                end
                
                obs_cls = bwconncomp(obs_tmp, conn);
                obs_labelmatrix = labelmatrix(obs_cls);
                obs_labelmatrix =  obs_labelmatrix(opt.spatial_def.inside);
                
                for i = 1:num_it
                    
                    r_tmp = zeros([opt.spatial_def.dim size(obs_map,2)]);
                    r_tmpstat = r_tmp;
                    for j = 1:size(obs_map,2)
                        r_tmp(opt.spatial_def.inside,j) = squeeze(r_map(i,:,j));
                        r_tmpstat(opt.spatial_def.inside,j) = squeeze(r_stat(i,:,j));
                    end
                    r_cls =  bwconncomp(r_tmp, conn);
                    if ~isempty(r_cls.PixelIdxList)  && strcmp(opt.clusterstatistic,'maxsize')
                        max_r_cls(i) = max(cellfun(@length,r_cls.PixelIdxList));
                    elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'maxsum')
                        tmp_cls = nan(1,length(r_cls.PixelIdxList));
                        for ii = 1:length(r_cls.PixelIdxList)
                            tmp_cls(ii) = sum(r_tmpstat(r_cls.PixelIdxList{ii}));
                        end
                        max_r_cls(i) = max(tmp_cls);
                    elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'max')
                        tmp_cls = nan(1,length(r_cls.PixelIdxList));
                        for ii = 1:length(r_cls.PixelIdxList)
                            tmp_cls(ii) = max(r_tmpstat(r_cls.PixelIdxList{ii}));
                        end
                        max_r_cls(i) = max(tmp_cls);
                    elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'wcm')
                        tmp_cls = nan(1,length(r_cls.PixelIdxList));
                        for ii = 1:length(r_cls.PixelIdxList)
                            tmp_cls(ii) = sum((r_tmpstat(r_cls.PixelIdxList{ii})-thresh).^1);
                        end
                        max_r_cls(i) = max(tmp_cls);
                    end
                end
                
                
            else
                
                %here we are simply looking for clusters w/o spatial  structure
                conn = conndef(ndims(obs_map),'max');
                obs_cls = bwconncomp(obs_map,conn);
                obs_labelmatrix = labelmatrix(obs_cls);
                
                for i = 1:num_it
                    
                    r_tmp = squeeze(r_map(i,:,:));
                    r_cls =  bwconncomp(r_tmp,conn);
                    r_tmpstat = squeeze(r_stat(i,:,:));
                    if ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'maxsize')
                        max_r_cls(i) = max(cellfun(@length,r_cls.PixelIdxList));
                    elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'maxsum')
                        tmp_cls = nan(1,length(r_cls.PixelIdxList));
                        for ii = 1:length(r_cls.PixelIdxList)
                            tmp_cls(ii) = sum(r_tmpstat(r_cls.PixelIdxList{ii}));
                        end
                        max_r_cls(i) = max(tmp_cls);
                   elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'max')
                        tmp_cls = nan(1,length(r_cls.PixelIdxList));
                        for ii = 1:length(r_cls.PixelIdxList)
                            tmp_cls(ii) = max(r_tmpstat(r_cls.PixelIdxList{ii}));
                        end
                        max_r_cls(i) = max(tmp_cls);
                    elseif ~isempty(r_cls.PixelIdxList) && strcmp(opt.clusterstatistic,'wcm')
                        tmp_cls = nan(1,length(r_cls.PixelIdxList));
                        for ii = 1:length(r_cls.PixelIdxList)
                            tmp_cls(ii) = sum((r_tmpstat(r_cls.PixelIdxList{ii})-thresh).^1);
                        end
                        max_r_cls(i) = max(tmp_cls);
                    end
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
                    obs_clstat(i) = sum(obs_stat(obs_cls.PixelIdxList{i}));
                    cluster_pvals(i) = ((sum(max_r_cls>=obs_clstat(i)))+1)/(num_it+1);
                end
            elseif strcmp(opt.clusterstatistic,'max')
                cluster_pvals = nan(1,length(obs_cls.PixelIdxList));
                obs_clstat = cluster_pvals;
                for i = 1:length(obs_cls.PixelIdxList)
                    obs_clstat(i) = max(obs_stat(obs_cls.PixelIdxList{i}));
                    cluster_pvals(i) = ((sum(max_r_cls>=obs_clstat(i)))+1)/(num_it+1);
                end
            elseif strcmp(opt.clusterstatistic,'wcm')
                cluster_pvals = nan(1,length(obs_cls.PixelIdxList));
                obs_clstat = cluster_pvals;
                for i = 1:length(obs_cls.PixelIdxList)
                    obs_clstat(i) = sum((obs_stat(obs_cls.PixelIdxList{i})-thresh).^1);
                    cluster_pvals(i) = ((sum(max_r_cls>=obs_clstat(i)))+1)/(num_it+1);
                end
                max_r_cls(i) = max(tmp_cls);
            end
            
            
            %save stuff
            stats.clusters = obs_cls.PixelIdxList;
            stats.clusterstat = obs_clstat;
            stats.clusterlabelmatrix = obs_labelmatrix;
            stats.clusterpvals = cluster_pvals;
            stats.randclusterstatmax = max_r_cls;
            stats.sigclusters = stats.clusters(stats.clusterpvals<alpha);
            stats.sigpvals = stats.clusterpvals(stats.clusterpvals<alpha);
            
            %save an indexing vector for time-resolved data
            if ismatrix(accuracy)
                stats.sigtime = zeros(1,size(accuracy,2));
                tpidx = cat(1,stats.clusters{stats.clusterpvals<alpha});
                stats.sigtime(tpidx) = 1;
            end
        end
        
        
    case {'mixed_of', 'mixed_fo'}
        
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