function [assemblies,triplets] = TolerantDefragmentSSA(origAssemblies,spikes,windowSize,threshold,nMin,verbose,tolerance)

if isempty(origAssemblies)
    assemblies = origAssemblies;
    return
end

if ~exist('tolerance','var') || isempty(tolerance), tolerance = 0; end

if verbose, tic; end
% origAssemblies = mergeAssemblies(origAssemblies')'; % This can be too time-consiming and should already have been performed
nUnits = size(origAssemblies,2);
sizes = sum(origAssemblies,2);
% start with the minimum sized assemblies
currentSize = min(sizes);
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
assemblies = cell2mat(cellData(~cellfun(@isempty,cellData)));

if verbose, display([datestr(clock) ': Merging new assemblies to avoid repetitions. toc=' addComma(round(toc))]); end
newAssemblies = unique(assemblies,'rows');
if verbose, display([datestr(clock) ': Done with preprocessing: starting with ' addComma(size(newAssemblies,1)) ' ' addComma(currentSize) '-member assemblies. toc=' addComma(round(toc))]); end
triplets = newAssemblies;
%%
change = true;
for currentSize = (currentSize):100 % It should technically be Infinity but 100 is used to prevent the function from running indefinitely in case of a problem
    assemblies = newAssemblies;
    setAside{currentSize,1} = logical(assemblies);
    if verbose
        display([datestr(clock) ': Finding ' addComma(currentSize+1) '-member cliques in ' addComma(size(assemblies,1)) ' ' addComma(currentSize) '-member assemblies. toc=' addComma(round(toc))]);
    end
    if size(assemblies,1)<2, break; end
    these = double(assemblies);
    
    overlap = @(x) x*x' + diag(nan(size(x,1),1));
    antioverlap = @(x) x*(1-x)' + diag(nan(size(x,1),1));
    nonExtendable = zeros(0,size(these,2));
    while change
        n0 = size(these,1);
        try
	  maxOverlap = max(overlap(these),[],2);
	  minAntioverlap = min(antioverlap(these),[],2);
	  nonExtendable = [nonExtendable; these(maxOverlap<2,:)];
	  these(maxOverlap<2,:) = [];
        end
        
        possibleMembers = sum(these,1)>=(currentSize-(tolerance>0));
        nonExtendable = [nonExtendable; these(sum(these(:,~possibleMembers),2)>0,:)];
        these(sum(these(:,~possibleMembers),2)>0,:) = [];
        
        if size(these,1)==n0 || isempty(these), change=false; end % quit loop
    end
    if isempty(these),break; end
    % Consider adding any of the possible members to any assembly
    nPossible = sum(possibleMembers);
    nAll = factorial(nPossible)/(factorial(nPossible-(currentSize+1))*factorial(currentSize+1));
    if verbose
        display([datestr(clock) ': Considering all ' addComma(min([nAll, size(these,1)*nPossible])) ' ways to extend the ' addComma(currentSize) '-member assemblies. toc=' addComma(round(toc))]);
    end
    % We need to do a loop for all possible cliques containing assembly (j)
    % (considering only assemblies after 'j')
    [cellCliques,cellHits] = deal(cell(size(these,1),1));
    
    try
        for j=1:size(these,1) % for assembly j
	  possibleExtensions = [ j+find(sum(these(j+1:end,these(j,:)>0),2)==currentSize-1)];
	  % remove neurons appearing too few times to possibly participate in a clique
	  possibleNeurons = sum(these([j;possibleExtensions],:),1)>=2;
	  possibleExtensions(any(these(possibleExtensions,~possibleNeurons),2)) = [];
	  possiblePartners = find(sum(these(:,possibleNeurons),2)==currentSize);
	  combinations = these(possibleExtensions,:); combinations(:,these(j,:)>0) = 1; combinations = unique(combinations,'rows');
	  
	  theseHits = cell(size(combinations,1),1); theseCliques = false(size(combinations,1),1);
	  for i=1:size(combinations,1)
	      ok = combinations(i,:)>0;
	      comb = SSA_SplitAssemblies(combinations(i,:),3);
	      if mean(~ismember(comb(:,ok),triplets(:,ok),'rows'))<=tolerance
		theseCliques(i) = true;
		hitID = sum(these(possibleExtensions,ok),2)==currentSize;
		theseHits{i} = possiblePartners(hitID);
	      end
	  end
	  cellCliques{j} = combinations(theseCliques,:); cellHits{j} = theseHits(theseCliques,:);
	  
	  if verbose && (rem(j,500)==0)
	      display([datestr(clock) ': So far found ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques after going through ' ...
		addComma(j) '/' addComma(size(these,1)) ' groups. toc=' addComma(round(toc))]);
	  end
	  
        end
    catch
        keyboard
    end
    %     end
    
    if verbose
        display([datestr(clock) ': Considering each of the ' addComma(sum(cellfun(@(x) size(x,1),cellCliques))) ' cliques if it satisfies the MUA criterion. toc=' addComma(round(toc))]);
    end    
    combinations = cell2mat(cellCliques);
    hits = cat(1,cellHits{:});
    if tolerance>0
        [~,u] = unique(combinations,'rows');
        combinations = combinations(u,:);
        hits = hits(u,:);
    end
    n = floor(size(combinations,1)/1000)*1000;
    pass = false(size(combinations));

    if threshold~=-Inf
        if verbose,tic; end
        for loop = 000:1000:n-1 % do this in batches so that results are periodically saved in the cell, rather than kept by each worker until the end
	  if n==0, continue; end
	  if verbose
	      display([datestr(clock) ': Loop ' addComma(loop) ' out of ' addComma(n) '. toc=' addComma(round(toc))]);
	  end
	  parfor j=(1:1000)+loop % length(cellCliques)
	      thisCombination = combinations(j,:);
	      pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
	      if verbose && rem(j,100)==0, display([datestr(clock) ': ' num2str(j)]); end
	  end
        end
        if verbose
	  display([datestr(clock) ': Last loop containing ' addComma(sum(cellfun(@(x) size(x,1),cellCliques((n+1):length(pass))))) ' cliques. toc=' addComma(round(toc))]);
        end
        verbose = verbose; % make this variable accessible to the workers
        parfor j=(n+1):size(combinations,1)% length(cellCliques)
	  thisCombination = combinations(j,:);
	  pass(j) = RateSSA2(thisCombination,spikes,windowSize,threshold,nMin);
	  if verbose && rem(j,100)==0, display([datestr(clock) ': cell passes ' num2str(j) ' computed.']); end
        end
    else % If the threshold is -Inf, don't check the first criterion
        pass = true(size(combinations,1),1);
    end
    
    these(cell2mat(hits(pass)),:) = [];
    newAssemblies = combinations(pass,:);
    if verbose,
        display([datestr(clock) ': ' addComma(sum(pass)) ' ' addComma(currentSize+1) '-member cliques found and saved. ' ...
	  addComma(size(assemblies,1)) ' ' addComma(currentSize) '-member assemblies did not participate in cliques (also saved). toc=' addComma(round(toc))]);
    end
    setAside{currentSize,1} = logical([these; nonExtendable]);
    setAside{currentSize,2} = logical(newAssemblies);
    date = datestr(clock); date(strfind(date,' ')) = '-';date(min(strfind(date,':'))) = 'h'; date(end-2:end) = [];
    
    if isempty(newAssemblies),break; end
    
    %     if verbose
    try save(['temp-defragmented-' date '-' num2str(currentSize) '.mat'],'setAside','-v7.3'); end
    try eval(['!rm ' 'temp-defragmented-' date0 '-' num2str(currentSize-1) '.mat']); end % remove last temp file
    date0 = date;
    %     end
end

try
    assemblies = cell2mat(setAside(~cellfun(@isempty,setAside(:,1)),1));
    assemblies = mergeAssemblies(assemblies')';
catch
    keyboard
end

try eval(['!rm ' 'temp-defragmented-' date0 '-' num2str(currentSize-1) '.mat']); end % remove temp file

function numOut = addComma(numIn) % for display purposes

% Found online, written by Ted Shultz. This function adds commas when displaying a number:
% EXAMPLE:
% numOut = addComma(1025620); % produces the string 1,025,620


jf=java.text.DecimalFormat; % comma for thousands, three decimal places
numOut= char(jf.format(numIn)); % omit "char" if you want a string out


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







