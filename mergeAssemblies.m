function assemblies = mergeAssemblies(assemblies, limit3)
% This function removes assemblies contained within others.
% 
% INPUT:
% - assemblies: M x P matrix with M the number of cells, and P the number
%   of assemblies.
%   --OR-- Nx1 cell, with N the different times scales. Each cell
%   contains an M x P matrix. 
% 
% OUTPUT:
% - mergedAssemblies: M x V matrix, with M the number of cells, and V the
%   number of unique matrices.
%

if iscell(assemblies)
    assemblies = cell2mat(assemblies');
end

assemblies=unique(assemblies','rows')';
commons = double(assemblies)' * double(assemblies); 
commons(1:length(commons)+1:end) = 0;
numNeurons = sum(assemblies);
isIncluded = bsxfun(@eq, commons,numNeurons);
assemblies(:,sum(isIncluded)>0) = [];
assemblies=unique(assemblies','rows')';

if ~exist('limit3','var')
    return
end
    
if limit3 == 1
    assemblies(:,sum(assemblies)<3) = [];
end

