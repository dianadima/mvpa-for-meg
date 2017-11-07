function [ data,labels ] = prepare_meg_data( data1, data2, varargin )
% Input: full path of files containing channel x time x trial data for each
%    condition. Data variable can have any name, but has to be the only variable 
%    in the file. Trial dimension should be the last one. Otherwise, specify the
%    trial dimension using a 3rd argument (i.e. 1, if data is trials*features).
% Output: classifier-ready concatenated data & label vector.

name = whos('-file', data1);
var = name(1).name;
load(data1, var);
eval(sprintf('data1 = %s;',var));

name = whos('-file',data2);
var = name(1).name;
load(data2, var);
eval(sprintf('data2 = %s;',var));

if ndims(data1) ~= ndims(data2)
    error('Both datasets need to have the same size and number of dimensions.')
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

