function [assemblies,triplets] = DefragmentSAC(origAssemblies,spikes,windowSize,threshold,nMin,verbose,tolerance)

% DefragmentSAC merges overlapping assemblies into larger candidate groups.
%
%   [assemblies,triplets] = DefragmentSAC(origAssemblies,spikes,...
%       windowSize,threshold,nMin,verbose,tolerance)
%
%   This function iteratively combines smaller overlapping assemblies into
%   larger assemblies when the merged groups satisfy co-activation
%   criteria. Larger assemblies are constructed from compatible smaller
%   subassemblies ("cliques").
%
%   REQUIRED INPUTS
%   origAssemblies - a binary matrix where each row corresponds to an
%                    assembly and each column corresponds to a neuron.
%
%   spikes         - a two-column [timestamp, unitID] matrix containing
%                    the list of spikes for each unit.
%
%   windowSize     - assembly timescale used to evaluate co-activation.
%
%   threshold      - statistical threshold used to evaluate candidate
%                    merged assemblies (in z-units)
%
%   nMin           - minimum number of assembly activations required for a
%                    merged assembly to be retained.
%
%   verbose        - set to true to display progress messages during
%                    execution.
%
%   tolerance      - (Default: 0) fraction of missing lower-order
%                    subassemblies allowed when evaluating candidate
%                    assemblies.
%
%   OUTPUT
%   assemblies     - binary matrix containing the final merged assemblies.
%
%   triplets       - binary matrix containing the initial lowest-order
%                    assembly fragments used during clique construction.
%
% Copyright (C) 2020-2026 by Ralitsa Todorova & Gabriel Makdah
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

if isempty(origAssemblies)
    assemblies = origAssemblies;
    return
end

if ~exist('tolerance','var') || isempty(tolerance), tolerance = 0; end
if ~exist('verbose','var') || isempty(verbose), verbose = false; end
if ~exist('threshold','var') || isempty(threshold), threshold = Inf; end
if ~exist('nMin','var') || isempty(nMin), nMin = 0; end

if verbose, tic; end

%% Prepare data for merging: split assemblies into fragments

% origAssemblies = mergeAssemblies(origAssemblies')'; % This can be too time-consiming and should already have been performed
nUnits = size(origAssemblies,2);
sizes = sum(origAssemblies,2);
% start with the minimum sized assemblies
currentSize = min(sizes);
% We need to break up bigger assemblies in smaller constituents:
maxSize = max(sizes);
assemblies = origAssemblies(sizes==currentSize,:);
if verbose, disp([datestr(clock) ': ' displayNumber(size(assemblies,1)) ' ' num2str(currentSize) '-member assemblies present. Splitting larger assemblies into smaller constituents. Tic.']); end
cellData{currentSize,1} = logical(assemblies);
for i=currentSize+1:maxSize
    largeAssemblies = origAssemblies(sizes==i,:);
    if verbose, disp([datestr(clock) ': Starting to split ' displayNumber(size(largeAssemblies,1)) ' ' displayNumber(i) '-member assemblies into groups of ' displayNumber(currentSize) ' members. toc=' displayNumber(round(toc))]); end
    currentData = cell(size(largeAssemblies,1),1);
    for k=1:size(largeAssemblies,1)
        variants = nchoosek(find(largeAssemblies(k,:)),currentSize);
        fragment = false(size(variants,1),nUnits);
        fragment(sub2ind(size(fragment),repmat((1:size(variants,1))',currentSize,1),variants(:))) = true;
        currentData{k} = fragment;
    end
    cellData{i,1} = cell2mat(currentData);
    if verbose && toc>1, disp([datestr(clock) ': ' displayNumber(i) '-member assemblies split into groups of ' displayNumber(currentSize) ' members. toc=' displayNumber(round(toc))]); end
end
assemblies = cell2mat(cellData(~cellfun(@isempty,cellData)));

if verbose, disp([datestr(clock) ': Merging new assemblies to avoid repetitions. toc=' displayNumber(round(toc))]); end
newAssemblies = unique(assemblies,'rows');
if verbose, disp([datestr(clock) ': Done with preprocessing: starting with ' displayNumber(size(newAssemblies,1)) ' ' displayNumber(currentSize) '-member assemblies. toc=' displayNumber(round(toc))]); end
triplets = newAssemblies;

%% Merge the currentSize-member assemblies into (currentSize+1)-member assemblies when the addmissible

for currentSize = (currentSize):100 % It should technically be Infinity but 100 is used to prevent the function from running indefinitely in case of a problem
    assemblies = newAssemblies; % the "new assemblies" of the previous loop are the assemblies to merge on this loop
    setAside{currentSize,1} = logical(assemblies);
    if verbose, disp([datestr(clock) ': Finding ' displayNumber(currentSize+1) '-member cliques in ' displayNumber(size(assemblies,1)) ' ' displayNumber(currentSize) '-member assemblies. toc=' displayNumber(round(toc))]); end
    if size(assemblies,1)<2, break; end % cannot merge a single assembly, we are done

    candidateAssembliesToMerge = double(assemblies);
    % Some assemblies are not even theoretically extendable (do not have enough partners with overlapping members).
    % Save time by skipping attempting to merge these assemblies:
    change = true;
    overlap = @(x) x*x' + diag(nan(size(x,1),1));
    antioverlap = @(x) x*(1-x)' + diag(nan(size(x,1),1));
    nonExtendable = zeros(0,size(candidateAssembliesToMerge,2));
    while change
        n0 = size(candidateAssembliesToMerge,1);
        try
            maxOverlap = max(overlap(candidateAssembliesToMerge),[],2);
            minAntioverlap = min(antioverlap(candidateAssembliesToMerge),[],2);
            nonExtendable = [nonExtendable; candidateAssembliesToMerge(maxOverlap<2,:)];
            candidateAssembliesToMerge(maxOverlap<2,:) = [];
        end
        possibleMembers = sum(candidateAssembliesToMerge,1)>=(currentSize-(tolerance>0));
        nonExtendable = [nonExtendable; candidateAssembliesToMerge(sum(candidateAssembliesToMerge(:,~possibleMembers),2)>0,:)];
        candidateAssembliesToMerge(sum(candidateAssembliesToMerge(:,~possibleMembers),2)>0,:) = [];
        % Repeat this step as taking assemblies away would reveal that other assemblies are also non extendable. Repeat until no new assemblies have
        % been takken out (i.e. no change in this loop)
        if size(candidateAssembliesToMerge,1)==n0 || isempty(candidateAssembliesToMerge), change=false; end % quit loop
    end
    if isempty(candidateAssembliesToMerge),break; end % If no candidate assemblies are left, we are done.

    % Consider adding any of the possible members to any assembly
    nPossible = sum(possibleMembers);
    if nPossible>currentSize, nAll = factorial(nPossible)/(factorial(nPossible-(currentSize+1))*factorial(currentSize+1)); else, nAll = 0; end
    if verbose, disp([datestr(clock) ': Considering all ' displayNumber(min([nAll, size(candidateAssembliesToMerge,1)*nPossible])) ' ways to extend the ' displayNumber(currentSize) '-member assemblies. toc=' displayNumber(round(toc))]); end

    % We will do a loop for all possible cliques containing assembly (j)
    % (considering only assemblies after 'j')
    [cellCliques,cellHits] = deal(cell(size(candidateAssembliesToMerge,1),1));
    % a "clique" would be a candidate merged assembly cellCliques{j} is a [nCandidateMergedAssemblies nUnits] matrix for cliques inluding assembly "j"
    % a "hit" would contain the indices of other assemblies in "candidateAssembliesToMerge" (in addition to "j") that can be merged together with "j" to form the respective clique

    for j=1:size(candidateAssembliesToMerge,1) % for assembly j
        possibleExtensions = [j+find(sum(candidateAssembliesToMerge(j+1:end,candidateAssembliesToMerge(j,:)>0),2)==currentSize-1)]; % possible partnering "candidateAssembliesToMerge"
        % remove neurons appearing too few times to possibly participate in a clique
        possibleNeurons = sum(candidateAssembliesToMerge([j;possibleExtensions],:),1)>=2;
        possibleExtensions(any(candidateAssembliesToMerge(possibleExtensions,~possibleNeurons),2)) = [];
        possiblePartners = find(sum(candidateAssembliesToMerge(:,possibleNeurons),2)==currentSize);
        combinations = candidateAssembliesToMerge(possibleExtensions,:); combinations(:,candidateAssembliesToMerge(j,:)>0) = 1; combinations = unique(combinations,'rows');
        theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
        for i=1:size(combinations,1)
            ok = combinations(i,:)>0;
            comb = SplitAssemblies(combinations(i,:),3);
            if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
                theseCliques(i) = true;
                hitID = sum(candidateAssembliesToMerge(possibleExtensions,ok),2)==currentSize;
                theseHits{i} = [j; possibleExtensions(hitID)];
            end
        end
        cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
        if verbose && (rem(j,500)==0)
            disp([datestr(clock) ': So far found ' displayNumber(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
                displayNumber(j) '/' displayNumber(size(candidateAssembliesToMerge,1)) ' groups. toc=' displayNumber(round(toc))]);
        end
    end

    if verbose, disp([datestr(clock) ': Considering each of the ' displayNumber(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques if it satisfies the MUA criterion. toc=' displayNumber(round(toc))]); end
    combinations = cell2mat(cellCliques);
    hits = cat(1,cellHits{:});
    if tolerance>0
        [~,u] = unique(combinations,'rows');
        combinations = combinations(u,:);
        hits = hits(u,:);
    end
    n = floor(size(combinations,1)/1000)*1000;
    pass = false(size(combinations,1),1);

    if threshold~=-Inf
        if verbose,tic; end
        for loop = 000:1000:n-1 % do this in batches so that results are periodically saved in the cell, rather than kept by each worker until the end
            if n==0, continue; end
            if verbose
                disp([datestr(clock) ': Loop ' displayNumber(loop) ' out of ' displayNumber(n) '. toc=' displayNumber(round(toc))]);
            end
            for j=(1:1000)+loop % length(cellCliques)
                thisCombination = combinations(j,:);
                pass(j) = RateSSA(thisCombination,spikes,windowSize,threshold,nMin);
                if verbose && rem(j,100)==0, disp([datestr(clock) ': ' num2str(j)]); end
            end
        end
        if verbose, disp([datestr(clock) ': Last loop containing ' displayNumber(size(combinations,1)-n) ' cliques. toc=' displayNumber(round(toc))]); end
        verbose = verbose; % make this variable accessible to the workers
        for j=(n+1):size(combinations,1)% length(cellCliques)
            thisCombination = combinations(j,:);
            pass(j) = RateSSA(thisCombination,spikes,windowSize,threshold,nMin);
            if verbose && rem(j,100)==0, disp([datestr(clock) ': cell passes ' num2str(j) ' computed.']); end
        end
    else % If the threshold is -Inf, don't check the first criterion
        pass = true(size(combinations,1),1);
    end

    candidateAssembliesToMerge(cell2mat(hits(pass)),:) = [];

    newAssemblies = combinations(pass,:);
    if verbose,
        disp([datestr(clock) ': ' displayNumber(sum(pass)) ' ' displayNumber(currentSize+1) '-member cliques found and saved. ' ...
            displayNumber(size(assemblies,1)) ' ' displayNumber(currentSize) '-member assemblies did not participate in cliques (also saved). toc=' displayNumber(round(toc))]);
    end
    setAside{currentSize,1} = logical([candidateAssembliesToMerge; nonExtendable]);
    setAside{currentSize,2} = logical(newAssemblies);
    date = datestr(clock); date(strfind(date,' ')) = '-';date(min(strfind(date,':'))) = 'h'; date(end-2:end) = [];

    if isempty(newAssemblies),break; end

    try save(['temp-defragmented-' date '-' num2str(currentSize) '.mat'],'setAside','-v7.3'); end
    try eval(['!rm ' 'temp-defragmented-' date0 '-' num2str(currentSize-1) '.mat']); end % remove last temp file
    date0 = date;
end

try
    assemblies = cell2mat(setAside(~cellfun(@isempty,setAside(:,1)),1));
    assemblies = mergeAssemblies(assemblies')';
catch
    keyboard
end

try eval(['!rm ' 'temp-defragmented-' date0 '-' num2str(currentSize-1) '.mat']); end % remove temp file


function numOut = displayNumber(numIn) % for display purposes
% Found online, written by Ted Shultz. This function adds commas when displaying a number:
% EXAMPLE:
% numOut = displayNumber(1025620); % produces the string 1,025,620
jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out


function [pass,zOthers] = RateSSA(assembly,spikes,windowSize,threshold,nMin)

% Only second criterion (ignore zItself)
% Second criterion is that each neuron should be active during assembly activations more than the general population
members = find(assembly);
zOthers = nan(length(members),1);
pass = true;
id = spikes(:,2);
nUnits = max(id);
if ~exist('nMin','var'), nMin = 1; end

for j=1:length(members) % no reason to test latest addition; we've already performed this step with it
    jSpikes = spikes(id==members(j));
    nSpikes = length(jSpikes);
    if nSpikes==0,  pass = false; continue; end
    nonj = 1:length(members); nonj(j) = [];
    otherMembers = members(nonj);
    ok = ismember(spikes(:,2),otherMembers);
    s = spikes(ok,:);

    code = zeros(nUnits,1); code(otherMembers) = 1:length(otherMembers); % renumber the second column
    s(:,2) = code(s(:,2));

    activity = commonIntervals(s,windowSize,length(otherMembers),length(otherMembers));
    activity = mean(activity,2);
    if sum(~isnan(activity))<nMin,  pass = false; continue; end
    activity = bsxfun(@plus,mean(activity,2),[-windowSize/2 windowSize/2]);

    %     nonmemberMUA = spikes(~ok,1);
    %     globalCount = sum(ExclusiveCountInIntervals(nonmemberMUA,activity));  % How many of the spikes are within the window distance around activations
    globalCount = sum(ExclusiveCountInIntervals(spikes(:,1),activity)) - sum(ExclusiveCountInIntervals(s(:,1),activity));  % How many of the spikes are within the window distance around activations

    count = sum(ExclusiveCountInIntervals(jSpikes,activity));
    zOthers(j) = zBinomialComparison(count,length(jSpikes),globalCount,sum(~ok));

    if any(zOthers(j)<threshold), pass = false; end
end







