function [pass,assembly,zItself,zOthers] = TestSAC(assembly,spikes,windowSize,threshold)

% This tests for same criteria, but we are confirming that all the members still pass in the newly extended assembly

% Initialise variables
members = find(assembly);
zItself = nan(length(members),length(members));
zOthers = nan(length(members),1);
pass = true;
% make sure "spikes" is a matrix of [timestamp id], and "spikeCell" is a cell of {timestamps1,timestamps2,...} for each neuron
if ismatrix(spikes), for i=1:max(spikes(:,2)), spikesCell{i} =  spikes(spikes(:,2)==i); end
elseif iscell(spikes), spikesCell = spikes; temp = spikesCell;
    for i=1:length(spikesCell), if size(spikesCell{i},2)>size(spikesCell{i},1), spikesCell{i} = spikeCell{i}'; end
        temp{i} = spikesCell{i}; temp{i}(:,2) = i;
    end
    spikes = sortrows(cat(1,temp{:}));
end
nUnits = length(spikesCell);

for j=1:length(members) % Consider each member neuron separately
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


