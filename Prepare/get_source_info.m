function [] = get_source_info( output_file, varargin)
%Gets indices of source clusters to go into a searchlight analysis.
%Inputs: output file name & path.
%Make sure to use same sourcemodel as for the beamformer analysis! 
%Optional inputs: 
% source_selection: selection of sources(all, inside or all90; default: inside)
% resolution: spatial resolution of searchlight(default: 10 mm);
% decoding_window: time window of interest for further analysis/plotting (default: [], whole time axis saved).
% sample_rate: of the dataset as it will be used in analysis (after any resampling); default 1200;
% sourcemodel_resolution: which FT template sourcemodel to use (default:10 mm). Can be 4,5,6,7.5,8,10. 

p = inputParser;
addParameter(p, 'source_selection', 'inside'); %all, inside or aal90
addParameter(p, 'resolution', 10); %searchlight resolution in mm
addParameter(p, 'decoding_window', [-0.1 0.9]); 
addParameter(p, 'sample_rate', 1200);
addParameter(p, 'sourcemodel_resolution',10);
parse(p, varargin{:});

[~, ftdir] = ft_version; %get FT directory

if p.Results.sourcemodel_resolution==7.5
    load(fullfile(ftdir, 'template', 'sourcemodel', 'standard_sourcemodel3d7point5mm'));
else
    load(fullfile(ftdir, 'template', 'sourcemodel', ['standard_sourcemodel3d' num2str(p.Results.sourcemodel_resolution) 'mm']));
end;
sourcemodel = ft_convert_units(sourcemodel, 'mm'); %#ok<NODEF>

if strcmp(p.Results.source_selection, 'inside')
    voxels = sourcemodel.pos(sourcemodel.inside,:); 
elseif strcmp(p.Results.source_selection,'all')
    voxels = sourcemodel.pos;
elseif strcmp(p.Results.source_selection, 'aal90')
    roi_idx = source_to_aal('standard_sourcemodel3d10mm');
    roi_idx = unique(cell2mat(roi_idx));
    voxels = sourcemodel.pos(sourcemodel.inside,:); 
    voxels = voxels(roi_idx,:);
else
    error('Source selection can be: all, inside or aal90.')
end;

res = p.Results.resolution; %searchlight resolution in mm
source_idx = cell(length(voxels),1);

for id = 1:length(voxels)
    source_idx{id} = find(voxels(:,1)>=voxels(id,1)-res&voxels(:,1)<=voxels(id,1)+res&voxels(:,2)>=voxels(id,2)-res&voxels(:,2)<=voxels(id,2)+res&voxels(:,3)>=voxels(id,3)-res&voxels(:,3)<=voxels(id,3)+res);
end;
fprintf('\nCalculated %d cluster indices...', length(source_idx));

time = p.Results.decoding_window(1):(1/p.Results.sample_rate):p.Results.decoding_window(2);
fprintf('\nTime axis goes from %f to %f...', time(1), time(end));

save(output_file, 'source_idx', 'time');

end

