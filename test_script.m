function [assemblies,allActivations,nMin] = SimSpikeAssemblies(spikes,windowSize,threshold,nMin,verbose,maxSize)

% same as CoActivityDetector
% "spikes" [sorted]
% windowSize

if isempty(spikes), assemblies = []; allActivations = []; nMin = 0; return; end
if ~exist('nMin','var') || isempty(nMin),nMin = max(10,(spikes(end,1)-spikes(1))/60/5); end % at least once every 5 minutes (and a minimum of 10 times)
if ~exist('threshold','var'),threshold = 5; end
if ~exist('verbose','var'),verbose = false; end
if ~exist('maxSize','var'),maxSize = 50; end% how many cycles to

%% Start with assembly sizes of 1

nUnits = max(spikes(:,2));
[final,assemblies] = deal(zeros(0,nUnits));
allActivations = {};

nUnits = max(spikes(:,2));

pairs = combnk(1:nUnits,2); % We start off with all possible pairs. The first cycle will consider all 3-cell combinations.
lArray = false(length(pairs),nUnits); % lArray stands for "logical array", where each line is an assembly and each column is a neuron, 
% so [1 1 0 1 0 0;...] means neurons 1, 2 and 4 together form the first assembly
lArray(sub2ind(size(lArray),(1:size(lArray,1))',pairs(:,1)))=true;
lArray(sub2ind(size(lArray),(1:size(lArray,1))',pairs(:,2)))=true;

% Cycle to add assemblies to them
if verbose, tic; end
% try
for cycle = 1:maxSize-2
	if isempty(lArray), break, end
	% Initialise variablesraly
	
	lArray0 = lArray;
	nPrevStep = size(lArray0,1);
	zArray = nan(nPrevStep,nUnits);
	if verbose
		'starting...'
	end
	
	% compute z-values (confidence that we should add a given cell to the assembly)
	for i=1:nPrevStep
	    toc0=toc;
	    members = find(lArray0(i,:));
	    zArray(i,:) = helper_function(spikes,members,windowSize,threshold,nMin,1);
	    ttt(i,cycle) = toc-toc0;
	end
	tocs(cycle,1) = toc;
	% add the cells with z passing the threshold
	[iID,neuronID] = find(zArray>threshold);
	nMaxAssemblies = sum(zArray(:)>threshold); % initialise array given the maximum possible size (size before merging)
	if nMaxAssemblies>0
		lArray = lArray0(iID,:); % take the corresponding assemblies of the old array
		lArray(sub2ind(size(lArray),(1:nMaxAssemblies)',neuronID))=true;
		ok = false(nMaxAssemblies,1);
		for j=1:nMaxAssemblies
			ok(j) = VerifycAssembly(lArray(j,:),spikes,windowSize,threshold,neuronID(j));
		end
		lArray(~ok,:) = [];
		noBranches = Accumulate(iID,ok,'size',nPrevStep)==0; % no new branches
	else
		noBranches = ones(nPrevStep,1);
		ok = 0;
		lArray = [];
	end
	if cycle>1
		% consider adding another member to the dead-branch assemblies (possibly adding 2 to assembly [1 30 31 35]) as only higher values have been considered for
		zArray = nan(nPrevStep,nUnits);
		for i=1:nPrevStep
		    toc0=toc;
			if noBranches(i)==0, continue; end % if no larger-id members can be added
			members = find(lArray0(i,:));
			zArray(i,:) = helper_function(spikes,members,windowSize,threshold,nMin,2); % consider smaller-id members again
			ttt2(i,cycle) = toc-toc0;
		end
		% add the cells with z passing the threshold
		[iID,neuronID] = find(zArray>threshold);
		nMaxAssemblies2 = sum(zArray(:)>threshold); % initialise array given the maximum possible size (size before merging)
		if nMaxAssemblies2>0
			lArray2 = lArray0(iID,:); % take the corresponding assemblies of the old array
			lArray2(sub2ind(size(lArray2),(1:nMaxAssemblies2)',neuronID(:)))=true;
			ok = false(nMaxAssemblies2,1);
			for j=1:nMaxAssemblies2
				ok(j) = VerifycAssembly(lArray2(j,:),spikes,windowSize,threshold,neuronID(j));
			end
			lArray2(~ok,:) = [];
			noBranches2 = Accumulate(iID,ok,'size',nPrevStep)==0; % no new branches
		else
			noBranches2 = ones(nPrevStep,1);
			ok = 0;
			lArray2 = [];
		end
		lArray = [lArray;lArray2];
		noBranches = noBranches & noBranches2;
		final = [final; lArray0(noBranches,:)];
	else
		lArray2 = [];
		nMaxAssemblies2 = 0;
	end
	lArray = unique(lArray,'rows','stable');
	if verbose
		[cycle size(lArray,1) size(lArray2,1) nMaxAssemblies nMaxAssemblies2 floor(toc)]
	end
	if cycle==maxSize-2 % if we have reached user-set limit of expansion
		final = [final; lArray]; % add the assemblies we would otherwise expend
	end
	tocs(cycle,2) = toc;
end

assemblies = final(sum(final,2)>2,:);
assemblies = unique(assemblies,'rows','stable');

% remove assemblies contained within bigger assemblies (e.g. 1 3 within 1 2 3)
origAssemblies = assemblies;
[assemblies,~,idx] = unique(assemblies,'rows');
commons = double(assemblies) * double(assemblies)' ;
commons(1:length(commons)+1:end) = 0;
numNeurons = sum(assemblies,2);
isIncluded = bsxfun(@eq, commons,numNeurons);
% Change the indices so they reflect the bigger assembly (which will be retained)
bad = (sum(isIncluded,2)>0);
assemblies(bad,:) = [];

nAssemblies = size(assemblies,1);
% catch
% 	keyboard
% end

%% Save activations

if nargout<2
	return
end
for i=1:size(assemblies,1)
	ok = ismember(spikes(:,2),find(assemblies(i,:)));
	code = cumsum(assemblies(i,:))';
	allActivations{i,1} = commonIntervals([spikes(ok,1) code(spikes(ok,2))], windowSize, sum(assemblies(i,:)), sum(assemblies(i,:)));
	%     allActivations{i,1}(:,2)=i;
end






















