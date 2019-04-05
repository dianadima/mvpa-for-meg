function [roi_idx, labels] = source_to_aal(sourcemodel)
%Outputs indices of inside sources in fieldtrip beamformer grid that fall within each of 90 AAL ROIs.
%Inputs: - sourcemodel: sourcemodel or cortical sheet or full path to it (e.g. 'fieldtrip/template/sourcemodel/standard_sourcemodel3d10mm').
%
%DC Dima 2017 (diana.c.dima@gmail.com)

[~, ftdir] = ft_version; %get FT directory
atlas = ft_read_atlas(fullfile(ftdir, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii')); %load AAL atlas
atlas = ft_convert_units(atlas, 'mm');

if ischar(sourcemodel)
    [~,~,ext] = fileparts(sourcemodel);
    if strcmp(ext, '.mat')
        load(sourcemodel); %load FT default sourcemodel
    else
        sourcemodel = ft_read_headshape(sourcemodel);
        sourcemodel.inside = true(size(sourcemodel.pos,1),1);
    end
end
sourcemodel = ft_convert_units(sourcemodel,'mm');

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodel2 = ft_sourceinterpolate(cfg, atlas, sourcemodel); %interpolate grid to atlas
inside = find(sourcemodel.inside);        %we want indices for the inside sources only (3,294 for a 10 mm grid)

roi_idx = cell(1,90);

for i = 1:90                                %loop through ROIs
    indx = find(sourcemodel2.tissue==i);
    roi_idx{i} = zeros(1,length(indx));
    
    for j = 1:length(indx)
        roi_idx{i}(j) = find(inside==indx(j)); %store indices
    end
end

labels = sourcemodel2.tissuelabel(1:90);

if isfield(sourcemodel, 'tri') %if it's a cortical surface, remove deep structures
    roi_idx(71:80) = [];
    labels(71:80) = [];
end

end
