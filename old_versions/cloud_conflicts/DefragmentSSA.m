function assemblies = DefragmentSSA(origAssemblies,spikes,windowSize,threshold,nMin,verbose)

if isempty(origAssemblies)
	assemblies = origAssemblies;
	return
end

if verbose, tic; end
origAssemblies = mergeAssemblies(origAssemblies')';
nUnits = size(origAssemblies,2);
sizes = sum(origAssemblies,2);
% start with the minimum sized assemblies
currentSize = min(sizes);
% We need to break up bigger assemblies in smaller constituents:
maxSize = max(sizes);
assemblies = origAssemblies(sizes==currentSize,:);
for i=currentSize+1:maxSize
	largeAssemblies = origAssemblies(sizes==i,:);
	for k=1:size(largeAssemblies,1)
		variants = nchoosek(find(largeAssemblies(k,:)),currentSize);
		fragment = false(size(variants,1),nUnits);
		fragment(sub2ind(size(fragment),repmat((1:size(variants,1))',currentSize,1),variants(:))) = true;
		assemblies = [assemblies; fragment];
	end
end

newAssemblies = mergeAssemblies(assemblies')';
%%
for currentSize = currentSize:maxSize*2
	assemblies = newAssemblies;
	if verbose
		[currentSize size(assemblies,1) floor(toc)]
	end
	setAside{currentSize,1} = assemblies;
	if size(assemblies,1)<2, break; end
	these = double(assemblies);
	change = true;
	
	overlap = @(x) x*x' + diag(nan(size(x,1),1));
	antioverlap = @(x) x*(1-x)' + diag(nan(size(x,1),1));
	while change
		n0 = size(these,1);
		maxOverlap = max(overlap(these),[],2);
		minAntioverlap = min(antioverlap(these),[],2);
		q = [maxOverlap minAntioverlap sum(bsxfun(@eq,overlap(these),maxOverlap) & bsxfun(@eq,antioverlap(these),minAntioverlap),2)];
		
		these(maxOverlap<2,:) = [];
		
		possibleMembers = sum(these,1)>=currentSize;
		these(sum(these(:,~possibleMembers),2)>0,:) = [];
		
		if size(these,1)==n0 || isempty(these), change=false; end % quit loop
	end
	if isempty(these),break; end
	% Consider adding any of the possible members to any assembly
	nPossible = sum(possibleMembers);
	nAll = factorial(nPossible)/(factorial(nPossible-(currentSize+1))*factorial(currentSize+1));
	if nAll<size(these,1)*nPossible
		indices = nchoosek(find(possibleMembers),currentSize+1);
		combinations = false(size(indices,1),nUnits);
		combinations(sub2ind(size(combinations),ndgrid(1:size(indices,1),1:size(indices,2)),indices)) = true;
	else
		combinations = repmat(these>0.5,sum(possibleMembers),1);
		combinations(sub2ind(size(combinations),(1:size(combinations,1))',repelem(find(possibleMembers)',size(these,1)))) = true;
		% remove any repeating combinations
		combinations = unique(combinations(sum(combinations,2)>currentSize,:),'rows');
	end
	
	% Consider each combination. Is it a clique?
	cliques = false(size(combinations,1),1); hits = cell(size(combinations(:,1)));
	parfor i=1:size(combinations,1)
		fragments = repmat(combinations(i,:),currentSize+1,1);
		fragments(sub2ind(size(fragments),(1:currentSize+1),find(combinations(i,:))))=0;
		% if every fragment is found in the detected assemblies: 
		% (used to be ~any(~ismember(fragments,these,'rows')) but its faster to demand a "currentSize" overlap between the assemblies
		if min(max(these*fragments'))==currentSize
			cliques(i) = true;
			hits{i,1} = find(ismember(assemblies,fragments,'rows'));
		end
	end
	
	pass = cliques;
	new = find(cliques ); new = new(max(combinations(cliques,:)*origAssemblies',[],2)~=currentSize+1); % otherwise assembly was already present in the original set
	newpass = false(length(new),1);
	parfor j=1:length(new)
		i = new(j);
		newpass(j) = RateSSA2(combinations(i,:),spikes,windowSize,threshold,nMin);
	end
	pass(new) = newpass;
	
	assemblies(unique(cell2mat(hits(find(pass)))),:) = [];
	newAssemblies = combinations(pass,:);
	
	setAside{currentSize,1} = logical(assemblies);
	setAside{currentSize,2} = newAssemblies;
	
	date = datestr(clock); date(strfind(date,' ')) = '-';date(min(strfind(date,':'))) = 'h'; date(end-2:end) = [];
	save(['/home/raly/projects/Marco/temp-' date '.mat'],'setAside','origAssemblies','newAssemblies');
	if isempty(newAssemblies),break; end
	
end

assemblies = cell2mat(setAside(~cellfun(@isempty,setAside(:,1)),1));
assemblies = mergeAssemblies(assemblies')';




function [pass,zOthers] = RateSSA2(assembly,spikes,windowSize,threshold,nMin)

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
    if any(activity(:)==0); activity(activity==0)=nan; activity = nanmean(activity,2); % in case we run into the bug Gabriel found
    else
        activity = mean(activity,2);
	end
	if sum(~isnan(activity))<nMin,  pass = false; continue; end
  
    nonmemberMUA = spikes(~ok,1);
    globalCount = length(fastPETH(nonmemberMUA,activity,[-1 1]*windowSize/2));
    
    count = length(fastPETH(jSpikes,activity,[-1 1]*windowSize/2));
    zOthers(j) = zBinomialComparison(count,length(jSpikes),globalCount,length(nonmemberMUA));
    
    if any(zOthers(j)<threshold), pass = false; end
end







