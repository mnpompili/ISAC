function [assemblies,nMin] = SimSpikeAssemblies(spikes,windowSize,threshold,nMin,verbose,maxSize,minSize)

% input "spikes" is a [timestamp neuronID] matrix (sorted in time)
% input "windowSize" is the desired timescale of an assembly, i.e. how close in time spikes have to be to be considered an assembly activation
% input "threshold" is the threshold for the two statistical criteria to qualify a neuron belonging to an assembly.
% The threshold is provided in z units (standard deviations). We recommend using z=2.57, which for a single-tailed hypothesis (we only test in one diretion), corresponds to p~=0.005
% input "nMin" sets the minimum number of activations for an assembly to be considered (assemblies with fewer activations are discarded).
% This is to reduce computation time and not compute statistical criterai for irrelevant assemblies.
% By default, "nMin" is once every 5 minutes, or a minimum of 10 for recordings shorter than 50 minutes.
% input "verbose" is a logical value (true/false, false by default), where if verbose==true, program will display progress on every cycle (using tic/toc functions)
% input "maxSize" is an optional parameter of the largest assembly to consider. This way the user has the option to, for example, compute triplets only
% input "minSize" is an optional parameter of the smallest assembly to consider. By default, cell pairs are not considered assemblies (minSize=3 neurons)

if isempty(spikes), assemblies = []; allActivations = []; nMin = 0; return; end % return 0 assemblies in absence of spikes
% Default parameter values:
if ~exist('nMin','var') || isempty(nMin),nMin = max(10,(spikes(end,1)-spikes(1))/60/5); end % at least once every 5 minutes (and a minimum of 10 times)
if ~exist('threshold','var'),threshold = 5; end
if ~exist('verbose','var'),verbose = false; end
if ~exist('maxSize','var'),maxSize = max(spikes(:,2)); end% how many cycles to perform
if ~exist('minSize','var'),minSize = 3; end% by default, cell pairs are not considered assembly (triplets at a minimum)


%% Start with assembly sizes of 1

nUnits = max(spikes(:,2));
[final,assemblies] = deal(zeros(0,nUnits)); % variables containing the final assemblies

nUnits = max(spikes(:,2));

pairs = combnk(1:nUnits,2); % We start off with all possible pairs. Starting from these pairs, the first cycle will consider all 3-cell combinations (adding a cell to each pair)
lArray = false(length(pairs),nUnits); % lArray stands for "logical array", where each line is an assembly and each column is a neuron,
% so [1 1 0 1 0 0;...] would mean that neurons 1, 2 and 4 together form the first assembly
lArray(sub2ind(size(lArray),(1:size(lArray,1))',pairs(:,1)))=true;
lArray(sub2ind(size(lArray),(1:size(lArray,1))',pairs(:,2)))=true;

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
        % e,.g. consider adding 10 to [1 9] but not 8 (8<9). This avoids duplicates and each triplet is only tested once
        parfor i=1:nPrevStep
            members = find(lArray0(i,:));
            skipCell{i,1} = 1:max(members);
        end
    else
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
        parfor i=1:nPrevStep
            these = all_neuronID(all_iID==i);
            skipCell{i,1} = these((duplicates(all_iID==i)));
        end
    end
    
    if verbose
        ['starting... toc:' num2str(toc)]
    end
    % compute z-values (confidence that we should add a given cell to the assembly)
    spikesCell; % rappel?
    parfor i=1:nPrevStep % for each of the possible assemblies to extend
        members = find(lArray0(i,:));
        zArray(i,:) = helper_zToExtendAssembly(spikesCell,spikes,members,windowSize,threshold,nMin,skipCell{i}); % the test is performed in a helper function to keep the code more readable
        % When helper_zToExtendAssembly is called this way (last input, the "mode", is 1), it will only consider adding higher members of the assembly.
        % For example, an assembly with members [1 2 5], the first considered neuron is neuron 6, because neurons 3 and 4 were already considered on the last step.
        % This is for computational speed and to not spend every cycle considering the same neurons over an over again.
        % These neurons would be considered again below (in the "if cycle>1" statement)
    end
    if verbose
        ['confirming assemblies... toc:' num2str(toc)]
    end
    % add the cells with z values passing the threshold
    [assemblyID,neuronID] = find(zArray>=threshold);
    nExtendedAssemblies = sum(zArray(:)>=threshold); % the number of extensions
    if nExtendedAssemblies>0 % If there is at least 1 new assembly found this round
        lArray = lArray0(assemblyID,:); % take the corresponding assemblies of the old array
        lArray(sub2ind(size(lArray),(1:nExtendedAssemblies)',neuronID(:)))=true; % To each one, add the new qualified member
        passes = false(nExtendedAssemblies,1); % initialise passing variable
        parfor j=1:nExtendedAssemblies
            % Verify that this new extended assembly still passes the stitistical criterion for *every* member (not just the newly added one)
            passes(j) = helper_confirmAssembly(lArray(j,:),spikesCell,spikes,windowSize,threshold,neuronID(j));
        end
        lArray(~passes,:) = []; % The assemblies that don't pass the criteria for all their members get discarded
        % Compute a logical vector (row, one for each of the assemblies to extend) showing us which of them did not get extended this round (the sum of "pass"-es is 0)
        noBranches = Accumulate(assemblyID,passes,'size',nPrevStep)==0; % these are the assemblies with no branches
    else
        noBranches = ones(nPrevStep,1); % If no assemblies were extended, then they are all dead branches
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

% %% Save activations
%
% if nargout<2 % compute only if second argument is called
%     return
% end
% for i=1:size(assemblies,1)
%     ok = ismember(spikes(:,2),find(assemblies(i,:)));
%     code = cumsum(assemblies(i,:))';
%     allActivations{i,1} = commonIntervals_fast([spikes(ok,1) code(spikes(ok,2))], windowSize, sum(assemblies(i,:)), sum(assemblies(i,:)));
% end



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
if any(activity(:)==0); activity(activity==0)=nan; activity = nanmean(activity,2); % in case we run into a bug of the function where the assembly activation starts at 0
else
    activity = mean(activity,2);
end
if size(activity,1)<nMin % after cycle 1, such assemblies won't have qualified, but this avoids unneeded computations on cycle 1
    return
end
duration = sum(diff(ConsolidateIntervalsFast([activity-windowSize/2 activity+windowSize/2]),[],2)); % total duration of windows around the assembly activations

% Find incomplete activations of the assembly-to-be-extended
activityIncomplete = cell(nMembers,1); durationCountIncomplete = nan(nMembers,1); % each cell contains activations without a single member (e.g. activityIncomplete{1} contains activations of all members but the 1st member)
for without = 1:nMembers
    if nMembers==2
        activityIncomplete{without} =s(s(:,2)~=without,1); % when the assembly is just a pair of neurons, the "activations" excluding one neuron are ust the spikes of the other neuron
    else
        sIncomplete = s(s(:,2)~=without,:); % exclude spikes emitted by the member to be excludeds
        sIncomplete(sIncomplete(:,2)>without,2) = s(s(:,2)>without,2)-1'; % renumber the other neurons
        activityIncomplete{without} = commonIntervals_fast(sIncomplete,windowSize,nMembers-1,nMembers-1); % find incomplete activations
        % This next line detects and fixes a possible bug in "commonIntervals_fast" where the first column is an erroneous 0
        if any(activityIncomplete{without}(:)==0); activityIncomplete{without}(activityIncomplete{without}==0)=nan; activityIncomplete{without} = nanmean(activityIncomplete{without},2);
        else activityIncomplete{without} = mean(activityIncomplete{without},2); end
    end
    % Compute duration of incomplete activations
    durationCountIncomplete(without) = sum(diff(ConsolidateIntervalsFast([activityIncomplete{without}-windowSize/2 activityIncomplete{without}+windowSize/2]),[],2));
end

zOthers = zeros(nUnits,1); % First criterion (for each unit, does it fire significantly more with the assembly than the multi-unit activity fires with the assembly)
zItself = zeros(nUnits,nMembers); % Second criterion (for each unit, does it fire significantly more with the assembly than it does with each of the incomplete versions of the assembly)

% Control for MUA fluctuations
nonmemberMUA = spikes(~ok,1);
% How many of the global MUA spikes participate in activations:
globalCount = sum(CountInIntervals(activity,[nonmemberMUA-windowSize/2 nonmemberMUA+windowSize/2])>0);  % How many of the spikes are within the window distance around activations

% Consider adding any unit except those already members and those explictly excluded (duplicates excluded to save computation time)
neuronsToConsider = 1:nUnits;
neuronsToConsider([members(:); skip(:)]) = [];

for j=neuronsToConsider
    s = spikesCell{j};
    nSpikes = length(s);
    if isempty(s), continue; end
    
    % Are there enough activations containing this spike?
    count2 = sum(CountInIntervals(s,[activity-windowSize/2 activity+windowSize/2])>0); % How many activations are within the window distance around spikes
    if count2<nMin, continue; end % If the total number of activations including this neuron would be too low to consider the extended assembly, don't add this neuron and abort further tests
    
    % How many of this neuron's spikes participate in activations:
    count = sum(CountInIntervals(activity,[s-windowSize/2 s+windowSize/2])>0); % How many of the spikes are within the window distance around activations
    
    % Compute first statistical criterion (zOthers):
    zOthers(j) = zBinomialComparison(count,nSpikes,globalCount,length(nonmemberMUA));
    % zOthers is the proportion of spikes within activation for this neuron significantly higher than the equivalent proportion for the global multiunit activity
    
    if zOthers(j)>threshold % If the first criterion passes, compute the second one (zItself) as well
        for without = 1:nMembers
            countIncomplete = sum(CountInIntervals(activityIncomplete{without},[s-windowSize/2 s+windowSize/2])>0);% how many activations does the neuron participate in
            countIncomplete = countIncomplete / durationCountIncomplete(without)*duration; % normalise for duration (transform countIncomplete to a value comparable with the complete "count")
            zItself(j,without) = zBinomialComparison(count,nSpikes,countIncomplete,nSpikes); % perform a z-test (equivalent to a chi-square test)
            % if zItself(j,without)>threshold, then neuron j has significantly more of its spikes participating in a complete assembly activation
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
    if any(activity(:)==0); activity(activity==0)=nan; activity = nanmean(activity,2); % in case we run into the bug Gabriel found
    else
        activity = mean(activity,2);
    end
    duration = sum(diff(ConsolidateIntervalsFast([activity-windowSize/2 activity+windowSize/2]),[],2));
    
    % Now test for the two statistical criteria, starting from the second (stricter) one (zItself)
    activityIncomplete = cell(length(otherMembers),1); durationCountIncomplete = nan(length(otherMembers),1); % Compute incomplete activations of the assembly (excluding an additional neuron beyond the one we are considering)
    for without = 1:length(otherMembers)
        if length(otherMembers)==2
            activityIncomplete{without} =s(s(:,2)~=without,1);
        else
            sIncomplete = s(s(:,2)~=without,:);
            sIncomplete(sIncomplete(:,2)>without,2) = s(s(:,2)>without,2)-1';
            activityIncomplete{without} = commonIntervals_fast(sIncomplete,windowSize,length(otherMembers)-1,length(otherMembers)-1);
            % This next line detects and fixes a possible bug in "commonIntervals_fast" where the first column is an erroneous 0
            if any(activityIncomplete{without}(:)==0); activityIncomplete{without}(activityIncomplete{without}==0)=nan; activityIncomplete{without} = nanmean(activityIncomplete{without},2);
            else activityIncomplete{without} = mean(activityIncomplete{without},2); end
        end
        durationCountIncomplete(without) = sum(diff(ConsolidateIntervalsFast([activityIncomplete{without}-windowSize/2 activityIncomplete{without}+windowSize/2]),[],2));
    end
    
    % How many of this neuron's spikes participate in activations:
    count = sum(CountInIntervals(activity,[jSpikes-windowSize/2 jSpikes+windowSize/2])>0); % How many of the spikes are within the window distance around activations
    
    % How many of the global MUA spikes participate in activations:
    nonmemberMUA = spikes(~ok,1);
    globalCount = sum(CountInIntervals(activity,[nonmemberMUA-windowSize/2 nonmemberMUA+windowSize/2])>0);  % How many of the spikes are within the window distance around activations
    
    % First criterion
    zOthers(j) = zBinomialComparison(count,length(jSpikes),globalCount,length(nonmemberMUA));
    % zOthers is the proportion of spikes within activation for this neuron significantly higher than the equivalent proportion for the global multiunit activity
    if zOthers(j)<threshold, pass = false; assembly(members(j))=0; return; end % if any of the members does not pass any more, abort any further tests
    
    % Second criterion
    for without = 1:length(otherMembers)
        countIncomplete = sum(CountInIntervals(activityIncomplete{without},[jSpikes-windowSize/2 jSpikes+windowSize/2])>0);
        countIncomplete = countIncomplete / durationCountIncomplete(without)*duration; % normalise for duration (transform countIncomplete to a value comparable with the complete "count")
        zItself(j,nonj(without)) = zBinomialComparison(count,nSpikes,countIncomplete,nSpikes);
        % if zItself(j,without)>threshold, then neuron j has significantly more of its spikes participating in a complete assembly activation
        % than an activation of the assembly without one member ("incomplete" activations)
    end
    if any(zItself(j,:)<threshold), pass = false; assembly(members(j))=0; return; end % if any of the members does not pass any more, abort any further tests
end






















