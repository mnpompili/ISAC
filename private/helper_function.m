function theseZs = helper_function(spikes,members,windowSize,threshold,nMin,mode)

% mode 1 is adding larger neurons (with id past the last member)
% most 2 is the opposite (neurons with id smaller than the last member)

ok = ismember(spikes(:,2),members);
s = spikes(ok,:);
nMembers = length(members);
nUnits = max(spikes(:,2));
theseZs = nan(1,nUnits);
interval = [-1 1]*windowSize/2;

code = zeros(nUnits,1); code(members) = 1:nMembers; % renumber the second column
s(:,2) = code(s(:,2));
activity = commonIntervals(s,windowSize,nMembers,nMembers);
if any(activity(:)==0); activity(activity==0)=nan; activity = nanmean(activity,2); % in case we run into the bug Gabriel found
else
	activity = mean(activity,2);
end
if size(activity,1)<nMin % after cycle 1, such assemblies won't have qualified, but this avoids unneeded computations on cycle 1
	return
end
duration = sum(diff(ConsolidateIntervalsFast([activity-windowSize/2 activity+windowSize/2]),[],2));
activityWithout = cell(nMembers,1); durationWithout = nan(nMembers,1);
for without = 1:nMembers
	if nMembers==2
		activityWithout{without} =s(s(:,2)~=without,1);
	else
		sWithout = s(s(:,2)~=without,:);
		sWithout(sWithout(:,2)>without,2) = s(s(:,2)>without,2)-1';
		activityWithout{without} = commonIntervals(sWithout,windowSize,nMembers-1,nMembers-1);
		if any(activityWithout{without}(:)==0); activityWithout{without}(activityWithout{without}==0)=nan; activityWithout{without} = nanmean(activityWithout{without},2); % in case we run into the bug Gabriel found
		else activityWithout{without} = mean(activityWithout{without},2); end
	end
	durationWithout(without) = sum(diff(ConsolidateIntervalsFast([activityWithout{without}-windowSize/2 activityWithout{without}+windowSize/2]),[],2));
end
zItself = zeros(nUnits,nMembers);
zOthers = zeros(nUnits,1);

% Control for MUA fluctuations
nonmemberMUA = spikes(~ok,1);
globalCount = length(fastPETH(nonmemberMUA,activity,interval));
if mode==1
	neuronsToConsider = max(members)+1:nUnits;% to avoid repetitions, only test combinations where the member size goes up
else
	neuronsToConsider = 1:max(members); neuronsToConsider(members) = [];
end

for j=neuronsToConsider
	s = spikes(spikes(:,2)==j);
	nSpikes = length(s);
	if isempty(s), continue; end
	% i.e. how many activations does the neuron participate in:
	n = sum(abs(activity(FindClosest(activity,s))-s)<windowSize/2); % for how many activations is the nearest spike (to the activation) less than windowSize/2 away
	if n<nMin, continue; end
	for without = 1:nMembers
		nWithout = sum(abs(activityWithout{without}(FindClosest(activityWithout{without},s))-s)<windowSize/2);
		nWithout = nWithout / durationWithout(without)*duration; % normalise
		zItself(j,without) = zBinomialComparison(n,nSpikes,nWithout,nSpikes);
		% if zItself(j,without)>1.96, then neuron j has significantly more of its spikes participating in a complete assembly activation
		% than an activation of the assembly without one member (member "without" )
	end
	if min(zItself(j,:))>threshold
		count = length(fastPETH(s,activity,interval));
		zOthers(j) = zBinomialComparison(count,length(s),globalCount,length(nonmemberMUA));
	end
end
theseZs = zOthers; % only computed if the first one passes
