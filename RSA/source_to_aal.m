function [roi_idx] = source_to_aal(sourcemodelfile)
%Outputs indices of inside sources in fieldtrip beamformer grid that fall within each of 90 AAL ROIs.
%Inputs: - sourcemodelfile is filename of grid used in FT beamforming, i.e. standard_sourcemodel10mm.
%        - savefile is filename to store cell array containing indices (roi_idx).

atlas = ft_read_atlas('/cubric/software/MEG/fieldtrip-20161011/template/atlas/aal/ROI_MNI_V4.nii'); %load AAL atlas
load(['/cubric/software/MEG/fieldtrip-20161011/template/sourcemodel/' sourcemodelfile], 'sourcemodel'); %load sourcemodel
sourcemodel = ft_convert_units(sourcemodel,'mm');

cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodel2 = ft_sourceinterpolate(cfg, atlas, sourcemodel); %interpolate grid to atlas
a = 1:length(sourcemodel.inside);
inside = a(sourcemodel.inside);        %we want indices for the inside sources only (3,294 for a 10 mm grid)

roi_idx = cell(1,90);

for i = 1:90                                %loop through ROIs
    indx = find(sourcemodel2.tissue==i);
    roi_idx{i} = zeros(1,length(indx));
    
    for j = 1:length(indx)
        roi_idx{i}(j) = find(inside==indx(j)); %store indices
    end;
end;

end
