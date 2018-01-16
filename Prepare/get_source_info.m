function [source_idx, labels] = get_source_info(varargin)
%Gets indices of source clusters to go into a searchlight analysis.
%Inputs: output file name & path.
%Make sure to use same sourcemodel as for the beamformer analysis! 
%Optional inputs: 
% source_selection: selection of sources(all, inside or aal90; default: inside)
% resolution: spatial resolution of searchlight(default: 10 mm);
% sourcemodel: which FT template sourcemodel to use (default:10 mm). Can be 4,5,6,7.5,8,10. 
[~, ftdir] = ft_version; %get FT directory

p = inputParser;
addParameter(p, 'source_selection', 'inside'); %all, inside or aal90
addParameter(p, 'resolution', 10); %searchlight resolution in mm
addParameter(p, 'sourcemodel', fullfile(ftdir, 'template', 'sourcemodel', 'standard_sourcemodel3d10mm.mat'));
parse(p, varargin{:});

if ischar(p.Results.sourcemodel)
    [~,~,ext] = fileparts(p.Results.sourcemodel);
    if isempty(ext) || strcmp(ext, '.mat')
        load(p.Results.sourcemodel);
    else
        sourcemodel = ft_read_headshape(p.Results.sourcemodel);
    end;
else
    sourcemodel = p.Results.sourcemodel;
end;

if ~isfield(sourcemodel, 'inside')
    sourcemodel.inside = true(size(sourcemodel.pos,1),1);
end;
sourcemodel = ft_convert_units(sourcemodel, 'mm');

if strcmp(p.Results.source_selection, 'aal90')
    
    [source_idx, labels] = source_to_aal(sourcemodel);
    
else
    if strcmp(p.Results.source_selection, 'inside')
        voxels = sourcemodel.pos(sourcemodel.inside,:);
    elseif strcmp(p.Results.source_selection,'all')
        voxels = sourcemodel.pos;
    else
        error('Source selection can be: all, inside or aal90.')
    end;
    
    res = p.Results.resolution; %searchlight resolution in mm
    source_idx = cell(length(voxels),1);
    
    for id = 1:length(voxels)
        source_idx{id} = find(voxels(:,1)>=voxels(id,1)-res&voxels(:,1)<=voxels(id,1)+res&voxels(:,2)>=voxels(id,2)-res&voxels(:,2)<=voxels(id,2)+res&voxels(:,3)>=voxels(id,3)-res&voxels(:,3)<=voxels(id,3)+res);
    end;
    fprintf('\nCalculated %d cluster indices...\n', length(source_idx));
    
end;


end

