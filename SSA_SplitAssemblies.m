function splitAssemblies = SSA_SplitAssemblies(origAssemblies,nMembers,verbose)


nUnits = size(origAssemblies,2);
sizes = sum(origAssemblies,2);
% start with the minimum sized assemblies
if ~exist('nMembers','var'), currentSize = min(sizes); else currentSize = nMembers; end
if ~exist('verbose','var'), verbose = false; end
% We need to break up bigger assemblies in smaller constituents:
maxSize = max(sizes);
assemblies = origAssemblies(sizes==currentSize,:);
if verbose, display([datestr(clock) ': ' addComma(size(assemblies,1)) ' ' currentSize '-member assemblies present. Splitting larger assemblies into smaller constituents. Tic.']); end
cellData{currentSize,1} = logical(assemblies);
for i=currentSize+1:maxSize
    largeAssemblies = origAssemblies(sizes==i,:);
    if verbose, display([datestr(clock) ': Starting to split ' addComma(size(largeAssemblies,1)) ' ' addComma(i) '-member assemblies into groups of ' addComma(currentSize) ' members. toc=' addComma(round(toc))]); end
    currentData = cell(size(largeAssemblies,1),1);
    for k=1:size(largeAssemblies,1)
        variants = nchoosek(find(largeAssemblies(k,:)),currentSize);
        fragment = false(size(variants,1),nUnits);
        fragment(sub2ind(size(fragment),repmat((1:size(variants,1))',currentSize,1),variants(:))) = true;
        currentData{k} = fragment;
    end
    cellData{i,1} = cell2mat(currentData);
    if verbose && toc>1, display([datestr(clock) ': ' addComma(i) '-member assemblies split into groups of ' addComma(currentSize) ' members. toc=' addComma(round(toc))]); end
end
splitAssemblies = cell2mat(cellData(~cellfun(@isempty,cellData)));
splitAssemblies = unique(splitAssemblies,'rows');
