%set toolbox paths
fprintf('\nThese scripts were tested with Matlab R2015a and Fieldtrip-20190205.\n')
toolbox_path = strrep(mfilename('fullpath'),'mvpa_setup','');

%add dependencies with subdirectories
addpath(genpath(fullfile(toolbox_path,'liblinear-2.11')));
addpath(genpath(fullfile(toolbox_path,'libsvm')));

%ensure only first-level directories (and not subfolders) are added to path
addpath(fullfile(toolbox_path, 'Prepare'));
addpath(fullfile(toolbox_path, 'Decode'));
addpath(fullfile(toolbox_path, 'Plot'));
addpath(fullfile(toolbox_path, 'Stats')); 
addpath(fullfile(toolbox_path, 'RSA'));

%subfolders contain less recommended versions of scripts - run them by
%navigating to the subfolders.