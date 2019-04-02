function [ data,labels ] = prepare_binary_data( data1, data2, varargin )
% Concatenate two datasets for a binary classification problem and create label vector.
%
% Input: data1, data2: Data variables, or full paths of .mat files containing data for each condition (e.g. output from read_meg_data). Data can be stored in matrices of any size.
%                                         Data variables in these files can have any name, but have to be the only variable in each file. 
%     OPTIONAL: trial dimension. Trial dimension is assumed to be the last one in the matrix (e.g. channel x time x trials), otherwise needs to be specified. 
%                                Example: 1 if data is trials x features.
%
% Output: data: classifier-ready concatenated data (feature1 x feature2 x .... x trials) 
%         labels: label vector.
%
% DC Dima 2017 (diana.c.dima@gmail.com)

if ischar(data1) && ischar(data2)
    
    name = whos('-file', data1);
    var = name(1).name;
    load(data1, var);
    eval(sprintf('data1 = %s;',var));
    
    name = whos('-file',data2);
    var = name(1).name;
    load(data2, var);
    eval(sprintf('data2 = %s;',var));
    
end

if isstruct(data1) && isstruct(data2) %Fieldtrip
    data1 = cat(3, data1.trial{:});
    data2 = cat(3, data2.trial{:});
end

if ndims(data1) ~= ndims(data2)
    error('Both datasets need to have the same number of dimensions.')
end;

if isempty(varargin) %trials are the last dimension
    trldim = ndims(data1); 
else
    trldim = varargin{1};
end;

labels = [ones(size(data1,trldim),1); 2*ones(size(data2,trldim),1)];
data = cat(trldim, data1, data2);

if trldim~=ndims(data)
    dim = 1:ndims(data); dim(length(dim)+1) = trldim; dim(trldim) = [];
    data = permute(data, dim); %trials are last dimension
end;

end

