function assemblies = ISAC(spikes, windowSize, varargin)

% ISAC finds groups of neurons that fire together (assemblies).
%
%   assemblies = ISAC(spikes, windowSize)
%   assemblies = ISAC(spikes, windowSize, 'Name', Value, ...)
%
%   This function finds statistically significant groups of neurons, called
%   "assemblies," that tend to fire spikes at the same time.
%
%   REQUIRED INPUTS
%   spikes        - a two-column [timestamp, unitID] matrix containing
%                   the list of  spikes for each unit
%   windowSize    - assembly timescale. This is the largest time gap
%                   allowed between spikes for them to be considered a part
%                   of the same assembly activation event.
%   <options>      optional list of property-value pairs (see table below)
%
%   =========================================================================
%      Properties       Values
%    -------------------------------------------------------------------------
%     'threshold'       (Default: 2.57) A statistical cutoff (in z-scores) to decide
%                       if a neuron belongs to an assembly. A higher value makes the
%                       detection more strict. 2.57 corresponds to a
%                       p-value of ~0.005. If two values are provided, the
%                       first would be used for the affinity criterion, and
%                       the second, for the unity criterion.
%     'nMin'            (Default: 1) The minimum number of times an assembly must be
%                       active to be included in the results.
%     'verbose'         (Default: false) Set to 'true' to display progress messages
%                       while the function is running.
%     'maxSize'         (Default: total number of neurons n) The largest number of
%                       neurons allowed in a single assembly.
%     'minSize'         (Default: 3) The smallest number of neurons allowed in a
%                       single assembly.
%     'groupID'         (Default: ones(n,1)) A vector specifying the
%                       a group ID for for each neuron. This can be a cell
%                       type, a recorded region, or any division of cells
%                       (e.g. reward-responsive, etc).
%     'constraints'     (Default: {}) A cell containing the group
%                       constraint for the desired assemblies. For example,
%                       {[2]} would indicate that the function should only
%                       look for assemblies containing at least one neurons
%                       from group 2 (see 'groupID' above), and {[1 2 3]}
%                       indicates that the function should look for
%                       cross-group assemblies containing members from
%                       groups 1, 2, and 3, together. Note you can provide
%                       multiple cells: {[1 2],4} will look for
%                       cross-group assemblies containing members from both
%                       groups 1 and 2, and also any assemblies containing
%                       members from group 4.
%    =========================================================================

%
%   OUTPUT
%   assemblies    - A structure that contains the detected neuron assemblies
%                   and the times they were active.
%
% Copyright (C) 2020-2026 by Ralitsa Todorova & Gabriel Makdah
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% If there are no spikes, return an empty result.
if isempty(spikes)
    assemblies = [];
    return;
end


% --- Input Parser: Manages function inputs and default values ---
p = inputParser;

% Define required inputs
addRequired(p, 'spikes', @(x) isnumeric(x) && size(x, 2) == 2);
addRequired(p, 'windowSize', @isscalar);

% Define optional name-value pairs
addParameter(p, 'threshold', sqrt(2) * erfcinv(0.01), @(x) isscalar(x) | isvector(x));
addParameter(p, 'nMin', 1, @isscalar);
addParameter(p, 'verbose', false, @islogical);
addParameter(p, 'maxSize', max(spikes(:,2)), @isscalar);
addParameter(p, 'minSize', 3, @isscalar);
addParameter(p, 'groupID', ones(max(spikes(:, 2)), 1), @(x) isvector(x) & length(x)==max(spikes(:,2)));
addParameter(p, 'constraints', {}, @iscell);

% Parse the inputs provided by the user
parse(p, spikes, windowSize, varargin{:});

% Assign the parsed inputs (or their defaults) to variables
threshold = p.Results.threshold;
nMin = p.Results.nMin;
verbose = p.Results.verbose;
maxSize = p.Results.maxSize;
minSize = p.Results.minSize;
groupID = p.Results.groupID;
constraints = p.Results.constraints;

if length(threshold)==1, threshold = ones(1,2) * threshold; end
%% Start with assembly sizes of 1

spikes = sortrows(spikes);
nUnits = max(spikes(:,2));
[final,assemblies] = deal(zeros(0,nUnits)); % variables containing the final assemblies

nUnits = max(spikes(:,2));

pairs = combnk(1:nUnits,2); % We start off with all possible pairs. Starting from these pairs, the first cycle will consider all 3-cell combinations (adding a cell to each pair)
% Initialize "lArray":
% "lArray" stands for "logical array",
% where each line is an assembly and each column is a neuron,
% so [1 1 0 1 0 0;...] would mean that neurons 1, 2 and 4 together form the first assembly
lArray = false(length(pairs),nUnits);
lArray(sub2ind(size(lArray),(1:size(lArray,1))',pairs(:,1)))=true;
lArray(sub2ind(size(lArray),(1:size(lArray,1))',pairs(:,2)))=true;

% Apply constraints (remove pairs that don't correspond to desired constraints)
if ~isempty(constraints)
    lArray = helper_ApplyConstraints(lArray,groupID,constraints);
end

%% Cycle to expand assemblies

if verbose, tic; end

% Pre-allocate the spikes of each neuron in a cell for optimised calling
for j=1:max(spikes(:,2))
    spikesCell{j,1} = spikes(spikes(:,2)==j);
end

for cycle = 1:maxSize-2 % 1st cycle is triplets, so the number of cycles should be (maxSize-2)
    if isempty(lArray), break, end % If no assemblies were found in the previous cycle, this terminates the algorithm (there is nothing more to extend)
    
    % Initialise variables
    lArray0 = lArray; % previously found assemblies are in "lArray0". They are to be extended into the new "lArray" of this cycle
    nPrevStep = size(lArray0,1); % How many assemblies were found in the previous step. Each one of these is a possible assembly to extend
    zArray = nan(nPrevStep,nUnits); % a matrix of z-values (statistical test), where zArray(i,j) is the z-value for extending the "i"-th assembly with neuron "j"
    if verbose
        'starting...' % report that the cycle is starting if "verbose" is true
    end
    if verbose
        ['precomputing... toc:' num2str(toc)]
    end
    
    if cycle==1
        % Pre-compute comparisons to skip
        % First cycle we test all triplets. We can only test triplets with members higher than the current members
        % e.g. consider adding 10 to [1 9] but not 8 (8<9). This avoids duplicates and each triplet is only tested once
        if ~isempty(constraints)
            skipCell = cell(nPrevStep,1);
        else
            parfor i=1:nPrevStep %PARFOR
                members = find(lArray0(i,:));
                skipCell{i,1} = 1:max(members);
            end
        end
    else
        try
            % (we don't need to test adding 1 to [2 3] (candidate assembly [1 2 3]) if we consider add 3 to [1 2] (candidate assembly also [1 2 3])
            [all_iID,all_neuronID] = meshgrid(1:nPrevStep,1:nUnits); all_iID = all_iID(:); all_neuronID = all_neuronID(:);
            all_lArray = lArray0(all_iID(:),:); % take the corresponding assemblies of the old array
            all_lArray(sub2ind(size(all_lArray),(1:size(all_lArray,1))',all_neuronID(:)))=true;
            notExtending = sum(all_lArray,2)==cycle+1; % only 2 members for first cycle, or 3 members for 2nd cycle and so on.
            % This is the assembly itself (adding neuron 1 to assembly [1 2] does not extend the assembly. Do no consider these combinations
            all_iID(notExtending) = []; all_neuronID(notExtending) = []; all_lArray(notExtending,:) = [];
            
            % Compute duplicates in a loop to avoid taking too much ram:
            done = false(size(all_lArray(:,1)));
            for j=1:size(all_lArray,2)
                ok = all_lArray(:,j) & ~done;
                [~,index,~] = unique(all_lArray(ok,:),'rows','stable');
                dd = true(sum(ok),1);
                dd(index) = false;
                duplicates(ok)=dd;
                done = done|ok;
            end
            skipCell = cell(nPrevStep,1);
            parfor i=1:nPrevStep %PARFOR
                these = all_neuronID(all_iID==i);
                skipCell{i,1} = these((duplicates(all_iID==i)));
            end
        catch % if the required ram causes an error
            skipCell = cell(nPrevStep,1);
            parfor i=2:nPrevStep %PARFOR
                l = lArray0(i,:);
                this = double(l);
                code = 1:nUnits; code(l) = [];
                overlap = sum(lArray0(1:i-1,l),2);
                these = sum(lArray0(overlap==cycle,~l),1)>0;
                skipCell{i,1} = code(these);
            end
        end
    end
    
    % Apply constraints (remove pairs that don't correspond to desired constraints)
    if ~isempty(constraints)
        skipCell = helper_ApplyConstraints(lArray,groupID,constraints,skipCell);
    end
    
    if verbose
        display(['starting... toc:' num2str(toc) 's.'])
    end
    % compute z-values (confidence that we should add a gixven cell to the assembly)
    spikesCell; % rappel?
    % PARFOR
    parfor i=1:nPrevStep % for each of the possible assemblies to extend
        members = find(lArray0(i,:));
        zArray(i,:) = helper_zToExtendAssembly(spikesCell,spikes,members,windowSize,threshold,nMin,skipCell{i}); % the test is performed in a helper function to keep the code more readable
        % When helper_zToExtendAssembly is called this way (last input, the "mode", is 1), it will only consider adding higher members of the assembly.
        % For example, an assembly with members [1 2 5], the first considered neuron is neuron 6, because neurons 3 and 4 were already considered on the last step.
        % This is for computational speed and to not spend every cycle considering the same neurons over an over again.
        % These neurons would be considered again below (in the "if cycle>1" statement)
    end
    % add the cells with z values passing the threshold
    [assemblyID,neuronID] = find(zArray>=threshold(2));
    nExtendedAssemblies = sum(zArray(:)>=threshold(2)); % the number of extensions
    if verbose
        display(['confirming ' num2str(nExtendedAssemblies) ' candidate assemblies... toc:' num2str(toc) 's.'])
    end
    if nExtendedAssemblies>0 % If there is at least 1 new assembly found this round
        lArray = lArray0(assemblyID,:); % take the corresponding assemblies of the old array
        lArray(sub2ind(size(lArray),(1:nExtendedAssemblies)',neuronID(:)))=true; % To each one, add the new qualified member
        passes = false(nExtendedAssemblies,1); % initialise passing variable
        parfor j=1:nExtendedAssemblies %PARFOR
            % Verify that this new extended assembly still passes the stitistical criterion for *every* member (not just the newly added one)
            passes(j) = helper_confirmAssembly(lArray(j,:),spikesCell,spikes,windowSize,threshold,neuronID(j));
        end
        lArray(~passes,:) = []; % The assemblies that don't pass the criteria for all their members get discarded
        % Compute a logical vector (row, one for each of the assemblies to extend) showing us which of them did not get extended this round (the sum of "pass"-es is 0)
        noBranches = Accumulate(assemblyID,passes,'size',nPrevStep)==0; % these are the assemblies with no branches
    else
        noBranches = true(nPrevStep,1); % If no assemblies were extended, then they are all dead branches
        lArray = [];
    end
    
    final = [final; lArray0(noBranches,:)]; % add the assemblies that cannot be further extended to the list of final assemblies
    
    lArray = unique(lArray,'rows','stable'); % remove any duplicates (e.g. from adding neuron 35 to assembly [1 2 30 31] in the "forward step" and also adding neuron 2 to assembly [1 30 31 35] in the "backward step")
    if verbose
        disp(sprintf(['cycle ' num2str(cycle) ' done after a total of ' num2str(toc) 's. Extended ' num2str(nExtendedAssemblies)  ' assemblies.'...
            ' \n Next cycle will try to extend ' num2str(size(lArray,1)) ' confirmed assemblies.']));
        %         [cycle size(lArray,1) size(lArray2,1) nExtendedAssemblies nExtendedAssemblies2 floor(toc)] % display progress
    end
    if cycle==maxSize-2 % if we have reached user-set limit of expansion
        final = [final; lArray]; % add the assemblies we would otherwise expand
    end
end

assemblies = final(sum(final,2)>=minSize,:); % only consider assemblies at or above the minumum size
assemblies = unique(assemblies,'rows','stable'); % discard duplicates

% remove assemblies contained within bigger assemblies (e.g. 1 3 within 1 2 3)
[assemblies,~,idx] = unique(assemblies,'rows');
commons = double(assemblies) * double(assemblies)' ;
commons(1:length(commons)+1:end) = 0;
numNeurons = sum(assemblies,2);
isIncluded = bsxfun(@eq, commons,numNeurons);
% Change the indices so they reflect the bigger assembly (which will be retained)
bad = (sum(isIncluded,2)>0);
assemblies(bad,:) = [];

%% === Helper functions: ===

function theseZs = helper_zToExtendAssembly(spikesCell,spikes,members,windowSize,threshold,nMin,skip)

% Initialise variables
ok = ismember(spikes(:,2),members);
s = spikes(ok,:);
nMembers = length(members);
nUnits = max(spikes(:,2));
theseZs = nan(1,nUnits);

% Find activations of the assembly-to-be-extended
code = zeros(nUnits,1); code(members) = 1:nMembers; % renumber the second column so it consists of members only
s(:,2) = code(s(:,2));
activity = commonIntervals_fast(s,windowSize,nMembers,nMembers);
activity = bsxfun(@plus,mean(activity,2),[-windowSize/2 windowSize/2]);

if size(activity,1)<nMin % after cycle 1, such assemblies won't have qualified, but this avoids unneeded computations on cycle 1
    return
end
overlap = diff(activity(:,1)); overlap(overlap>windowSize)=[]; correction = sum(overlap-windowSize);
duration  = windowSize*size(activity,1)+correction;

% Find incomplete activations of the assembly-to-be-extended
activityIncomplete = cell(nMembers,1); durationCountIncomplete = nan(nMembers,1); % each cell contains activations without a single member (e.g. activityIncomplete{1} contains activations of all members but the 1st member)
for without = 1:nMembers
    if nMembers==2
        activityIncomplete{without} = bsxfun(@plus,s(s(:,2)~=without,1),[-windowSize/2 windowSize/2]); % when the assembly is just a pair of neurons, the "activations" excluding one neuron are ust the spikes of the other neuron
    else
        sIncomplete = s(s(:,2)~=without,:); % exclude spikes emitted by the member to be excludeds
        sIncomplete(sIncomplete(:,2)>without,2) = s(s(:,2)>without,2)-1'; % renumber the other neurons
        activityIncomplete{without} = commonIntervals_fast(sIncomplete,windowSize,nMembers-1,nMembers-1); % find incomplete activations
        % This next line detects and fixes a possible bug in "commonIntervals_fast" where the first column is an erroneous 0
        if any(activityIncomplete{without}(:)==0); activityIncomplete{without}(activityIncomplete{without}==0)=nan; end
        activityIncomplete{without} = bsxfun(@plus,mean(activityIncomplete{without},2),[-windowSize/2 windowSize/2]);
    end
    % Compute duration of incomplete activations
    overlap = diff(activityIncomplete{without}(:,1)); overlap(overlap>windowSize)=[]; correction = sum(overlap-windowSize);
    durationCountIncomplete(without)  = windowSize*size(activityIncomplete{without},1)+correction;
end

zOthers = nan(nUnits,1); % First criterion (for each unit, does it fire significantly more with the assembly than the multi-unit activity fires with the assembly)
zItself = nan(nUnits,nMembers); % Second criterion (for each unit, does it fire significantly more with the assembly than it does with each of the incomplete versions of the assembly)

% Control for MUA fluctuations
% New code: it's actially much faster to compute the total spikes and subtract member ones than to set a new "nonmemberMUA" variable:
globalCount = sum(ExclusiveCountInIntervals(spikes(:,1),activity))-sum(ExclusiveCountInIntervals(s,activity));  % How many of the spikes are within the window distance around activations
% Old code:
% nonmemberMUA = spikes(~ok,1);
% % How many of the global MUA spikes participate in activations:
% globalCount = sum(ExclusiveCountInIntervals(nonmemberMUA,activity));  % How many of the spikes are within the window distance around activations

% Consider adding any unit except those already members and those explictly excluded (duplicates excluded to save computation time)
neuronsToConsider = 1:nUnits;
neuronsToConsider([members(:); skip(:)]) = [];

for j=neuronsToConsider
    s = spikesCell{j};
    nSpikes = length(s);
    if isempty(s), continue; end
    
    %     % Are there enough activations containing this spike?
    %     count2 = sum(ExclusiveCountInIntervals(s,activity)>0); % How many activations are within the window distance around spikes
    %     if count2<nMin, continue; end % If the total number of activations including this neuron would be too low to consider the extended assembly, don't add this neuron and abort further tests
    
    % How many of this neuron's spikes participate in activations:
    count = sum(ExclusiveCountInIntervals(s,activity)); % How many of the spikes are within the window distance around activations
    
    % Compute first statistical criterion (zOthers):
    zOthers(j) = zBinomialComparison(count,nSpikes,globalCount,sum(~ok));
    % zOthers is the proportion of spikes within activation for this neuron significantly higher than the equivalent proportion for the global multiunit activity
    
    if zOthers(j)>threshold(1) % If the first criterion passes, compute the second one (zItself) as well
        for without = 1:nMembers
            countIncomplete = sum(ExclusiveCountInIntervals(s,activityIncomplete{without}));% how many activations does the neuron participate in
            countIncomplete = countIncomplete / durationCountIncomplete(without)*duration; % normalise for duration (transform countIncomplete to a value comparable with the complete "count")
            zItself(j,without) = zBinomialComparison(count,nSpikes,countIncomplete,nSpikes); % perform a z-test (equivalent to a chi-square test)
            % if zItself(j,without)>threshold(2), then neuron j has significantly more of its spikes participating in a complete assembly activation
            % than an activation of the assembly without one member ("incomplete" activations)
        end
    end
end
theseZs = min(zItself,[],2); % If zItself exists, zOthers must have passed


function [pass,assembly,zItself,zOthers] = helper_confirmAssembly(assembly,spikesCell,spikes,windowSize,threshold,skip)

% This tests for same criteria, but we are confirming that all the members still pass in the newly extended assembly

% Initialise variables
members = find(assembly);
zItself = nan(length(members),length(members));
zOthers = nan(length(members),1);
pass = true;
nUnits = max(spikes(:,2));

for j=1:length(members) % Consider each member neuron separately
    if ismember(members(j),skip), continue; end % We skip the newly added neuron, because it was already confirmed to pass both criteria
    jSpikes = spikesCell{members(j)};
    nSpikes = length(jSpikes);
    
    % Find the activations of the assembly excluding the neuron we are considering:
    nonj = 1:length(members); nonj(j) = []; % indices of other neurons
    otherMembers = members(nonj);
    ok = ismember(spikes(:,2),otherMembers);
    s = spikes(ok,:); % "s" contains all spikes of the other members
    code = zeros(nUnits,1); code(otherMembers) = 1:length(otherMembers); % renumber the second column
    s(:,2) = code(s(:,2));
    activity = commonIntervals_fast(s,windowSize,length(otherMembers),length(otherMembers));
    activity = bsxfun(@plus,mean(activity,2),[-windowSize/2 windowSize/2]);
    overlap = diff(activity(:,1)); overlap(overlap>windowSize)=[]; correction = sum(overlap-windowSize);
    duration  = windowSize*size(activity,1)+correction;
    
    % Now test for the two statistical criteria, starting from the second (stricter) one (zItself)
    activityIncomplete = cell(length(otherMembers),1); durationCountIncomplete = nan(length(otherMembers),1); % Compute incomplete activations of the assembly (excluding an additional neuron beyond the one we are considering)
    for without = 1:length(otherMembers)
        if length(otherMembers)==2
            activityIncomplete{without} = bsxfun(@plus,s(s(:,2)~=without,1),[-windowSize/2 windowSize/2]); % when the assembly is just a pair of neurons, the "activations" excluding one neuron are ust the spikes of the other neuron
        else
            sIncomplete = s(s(:,2)~=without,:);
            sIncomplete(sIncomplete(:,2)>without,2) = s(s(:,2)>without,2)-1';
            activityIncomplete{without} = commonIntervals_fast(sIncomplete,windowSize,length(otherMembers)-1,length(otherMembers)-1);
            % This next line detects and fixes a possible bug in "commonIntervals_fast" where the first column is an erroneous 0
            if any(activityIncomplete{without}(:)==0); activityIncomplete{without}(activityIncomplete{without}==0)=nan; end
            activityIncomplete{without} = bsxfun(@plus,mean(activityIncomplete{without},2),[-windowSize/2 windowSize/2]);
        end
        overlap = diff(activityIncomplete{without}(:,1)); overlap(overlap>windowSize)=[]; correction = sum(overlap-windowSize);
        durationCountIncomplete(without)  = windowSize*size(activityIncomplete{without},1)+correction;
    end
    
    % How many of this neuron's spikes participate in activations:
    count = sum(ExclusiveCountInIntervals(jSpikes,activity)); % How many of the spikes are within the window distance around activations
    
    % How many of the global MUA spikes participate in activations:
    globalCount = sum(ExclusiveCountInIntervals(spikes(:,1),activity))-sum(ExclusiveCountInIntervals(s,activity));  % How many of the nonmember spikes are within the window distance around activations  % How many of the spikes are within the window distance around activations
    
    % First criterion
    zOthers(j) = zBinomialComparison(count,length(jSpikes),globalCount,sum(~ok));
    % zOthers is the proportion of spikes within activation for this neuron significantly higher than the equivalent proportion for the global multiunit activity
    if zOthers(j)<threshold(1), pass = false; assembly(members(j))=0; return; end % if any of the members does not pass any more, abort any further tests
    
    % Second criterion
    for without = 1:length(otherMembers)
        countIncomplete = sum(ExclusiveCountInIntervals(jSpikes,activityIncomplete{without}));% how many activations does the neuron participate in
        countIncomplete = countIncomplete / durationCountIncomplete(without)*duration; % normalise for duration (transform countIncomplete to a value comparable with the complete "count")
        zItself(j,nonj(without)) = zBinomialComparison(count,nSpikes,countIncomplete,nSpikes);
        % if zItself(j,without)>threshold, then neuron j has significantly more of its spikes participating in a complete assembly activation
        % than an activation of the assembly without one member ("incomplete" activations)
    end
    if any(zItself(j,:)<threshold(2)), pass = false; assembly(members(j))=0; return; end % if any of the members does not pass any more, abort any further tests
end

function varargout = helper_ApplyConstraints(assemblies,groupID,constraints,skipCell)
% Apply constraints (e.g. "SGs should be cross-structural").
% This will remove assemblies that do not follow the constraints.



grouped = nan(size(assemblies,1),max(groupID));
for i=1:max(groupID)
    grouped(:,i) = sum(assemblies(:,groupID==i),2);
end


% This function can be called either to confirm existing SGs,
% or to flag cells to skip when considering extending the assembly.
% Behavior in each case is different:

if nargin<4 % apply constraints to existing SGs
    pass = false(size(assemblies,1),1);
    for i=1:length(constraints)
        constraint = constraints{i}(:);
        constraint = Accumulate(constraint(:),1,'size',max(groupID));
        constraint = constraint(:)';
        missing = constraint-grouped; missing(missing<0) = 0;
        maxMissingAllowed = max([sum(constraint)-sum(assemblies(1,:)),0]);
        pass(sum(missing,2)<=maxMissingAllowed) = true;
    end
    filtered = assemblies(pass,:);
    varargout = {filtered,pass};
else
    missingCell = cell(1,length(constraints));
    for i=1:length(constraints)
        constraint = constraints{i}(:);
        constraint = Accumulate(constraint(:),1,'size',max(groupID));
        constraint = constraint(:)';
        missing = constraint-grouped; missing(missing<0) = 0;
        missingCell{i} = missing;
    end
    alreadyFullfulled = sum(cell2mat(cellfun(@(x) ~any(x,2),missingCell,'UniformOutput',false)),2)>0;
    avoidGroup = sum(cat(3,missingCell{:}),3)==0;
    avoidGroup(alreadyFullfulled,:) = false;
    groups = cell(1,max(groupID));
    for i=1:max(groupID), groups{i} = find(groupID(:)'==i); end
    
    rowKeep  = num2cell(avoidGroup, 2); % each row → 1×M logical mask
    result = cellfun(@(sk, rk) unique([groups{rk}, sk], 'sorted'), skipCell, rowKeep, 'UniformOutput', false);
    
    varargout = {result};
end





















