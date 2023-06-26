function [pass,assembly,zItself,zOthers] = VerifycAssembly(assembly,spikes,windowSize,threshold,skip)

members = find(assembly);
zItself = nan(length(members),length(members));
zOthers = nan(length(members),1);
pass = true;
nUnits = max(spikes(:,2));

for j=1:length(members)
	if ismember(j,skip), continue; end
	jSpikes = spikes(spikes(:,2)==members(j));
	nSpikes = length(jSpikes);
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
	duration = sum(diff(ConsolidateIntervalsFast([activity-windowSize/2 activity+windowSize/2]),[],2));
	activityWithout = cell(length(otherMembers),1); durationWithout = nan(length(otherMembers),1);
	for without = 1:length(otherMembers)
		if length(otherMembers)==2
			activityWithout{without} =s(s(:,2)~=without,1);
		else
			sWithout = s(s(:,2)~=without,:);
			sWithout(sWithout(:,2)>without,2) = s(s(:,2)>without,2)-1';
			activityWithout{without} = commonIntervals(sWithout,windowSize,length(otherMembers)-1,length(otherMembers)-1);
			if any(activityWithout{without}(:)==0); activityWithout{without}(activityWithout{without}==0)=nan; activityWithout{without} = nanmean(activityWithout{without},2); % in case we run into the bug Gabriel found
			else activityWithout{without} = mean(activityWithout{without},2); end
		end
		durationWithout(without) = sum(diff(ConsolidateIntervalsFast([activityWithout{without}-windowSize/2 activityWithout{without}+windowSize/2]),[],2));
	end
	
	% i.e. how many activations does the neuron participate in:
	n = sum(abs(activity(FindClosest(activity,jSpikes))-jSpikes)<windowSize/2); % for how many activations is the nearest spike (to the activation) less than windowSize/2 away
	for without = 1:length(otherMembers)
		nWithout = sum(abs(activityWithout{without}(FindClosest(activityWithout{without},jSpikes))-jSpikes)<windowSize/2);
		nWithout = nWithout / durationWithout(without)*duration; % normalise
		zItself(j,nonj(without)) = zBinomialComparison(n,nSpikes,nWithout,nSpikes);
		% if zItself(j,without)>1.96, then neuron j is significantly more likely to participate in a complete assembly activation
		% than an activation of the assembly without one member (member "without" )
	end
	if any(zItself(j,:)<threshold), pass = false; assembly(members(j))=0; return; end
	
	nonmemberMUA = spikes(~ok,1);
	globalCount = length(fastPETH(nonmemberMUA,activity,[-1 1]*windowSize/2));
	
	count = length(fastPETH(jSpikes,activity,[-1 1]*windowSize/2));
	zOthers(j) = zBinomialComparison(count,length(jSpikes),globalCount,length(nonmemberMUA));
	
	if any(zOthers(j)<threshold), pass = false; assembly(members(j))=0; return; end
end

