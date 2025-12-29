% Generate random spike times assuming poisson point process

% %% User defined variables
% sessionStart = 0;
% sessionEnd = 7200; % Session length (Seconds)
% firingRate_mean = 3; % Average firing rate of the cells
% firingRate_jitter = 4; % SD of the first run to jitter cell firing rate
% firingRate_jitter2 = 1; % Cells with negative firing rate are ran a second time with this SD
% numCells = 50; % Number of cells to simulate
% assembly_meanJitter = 0.03; % Average jitter (in seconds) to add to assembly activations
% assembly_stdJitter = 0; % SD of jitter to add to assembly activations
% assembly_minJitter = 0; % Minimum acceptable jitter (in seconds)
% assemblyMeanSpikes = 2; % Average number of spikes within an assemmbly activation
% sizeAssemblies = [3 4 5 6 7]; % Sizes of assemblies
% numAssemblies = [100 25 10 10 10]; % Number of assemblies, order the same as sizeAssemblies
% meanAActs = 100; % Average number of times assemblies are active
% stdAActs = 0; % SD of the first run to compute the number of times assemblies are active
% stdAActs2 = 0; % SD of second run if cells have negative firing rates
% generateAssemblyIdentity = 1; % 1: generate assemblies || 0: if you supply groundTruthAssemblies
% randomizeAssemblyMeanSpikes = 1; %  0: exactly assemblyMeanSpikes every assembly activation || 1: poisson point process with lambda = assemblyMeanSpikes || 2: poisson point process with lambda = cellFiringRate/ 2
% globalActivityVariation = 1; % Set to 1 to include global variations in cell firing rates
% setMeanCellActivity = 2; % set to 0 if you're supplying cell firing rates manually (in cells)
% generateSpikes = 1; % set to 0 if you're supplying spikes manually
% negativeCorrelations = 0; % set to 1 to include random negative correlations within the dataset

% sessionStart = 0;
% sessionEnd = 600; % Session length (Seconds)
% firingRate_mean = 3; % Average firing rate of the cells
% firingRate_jitter = 0; % SD of the first run to jitter cell firing rate
% firingRate_jitter2 = 0; % Cells with negative firing rate are ran a second time with this SD
% numCells = 10; % Number of cells to simulate
% assembly_meanJitter = 0.03; % Average jitter (in seconds) to add to assembly activations
% assembly_stdJitter = 0; % SD of jitter to add to assembly activations
% assembly_minJitter = 0; % Minimum acceptable jitter (in seconds)
% assemblyMeanSpikes = 2; % Average number of spikes within an assemmbly activation
% sizeAssemblies = [3 4]; % Sizes of assemblies
% numAssemblies = [4 1]; % Number of assemblies, order the same as sizeAssemblies
% meanAActs = 100; % Average number of times assemblies are active
% stdAActs = 0; % SD of the first run to compute the number of times assemblies are active
% stdAActs2 = 0; % SD of second run if cells have negative firing rates
% generateAssemblyIdentity = 1; % 1: generate assemblies || 0: if you supply groundTruthAssemblies
% randomizeAssemblyMeanSpikes = 2; %  0: exactly assemblyMeanSpikes every assembly activation || 1: poisson point process with lambda = assemblyMeanSpikes || 2: poisson point process with lambda = cellFiringRate/ 2
% globalActivityVariation = 0; % Set to 1 to include global variations in cell firing rates
% setMeanCellActivity = 2; % set to 0 if you're supplying cell firing rates manually (in cells)
% generateSpikes = 1; % set to 0 if you're supplying spikes manually
% negativeCorrelations = 0; % set to 1 to include random negative correlations within the dataset
% assemblyActivationTimestamps = {};

%% Add in assembly activations
numA = sum(numAssemblies);

if setMeanCellActivity == 0 || generateSpikes == 0
    if generateSpikes == 0
        spikes = Restrict(spikes, [sessionStart sessionEnd]);
        cells = spikesToCells(spikes);
    end
    numCells = length(cells);
    lambda = cellfun(@length,cells(1:numCells)) / (sessionEnd-sessionStart);
elseif setMeanCellActivity == 1
    lambda = randn(numCells,1) * firingRate_jitter + (firingRate_mean);
elseif setMeanCellActivity == 2
    lambda = randn(numCells,1) * firingRate_jitter + (firingRate_mean);
    while any(lambda < 0.5)
        lambda(lambda < 0.5) = randn(sum(lambda < 0.5),1) * firingRate_jitter2 + (firingRate_mean);
    end
end
% 

numPoints = floor(lambda * sessionEnd);
indices = [0; cumsum(numPoints)];

if generateSpikes == 1
    spikes = [];
    for i = 1:numCells
        tempSpikes = cumsum(-log(rand(numPoints(i),1))/lambda(i), 1);
        tempSpikes = tempSpikes./max(tempSpikes).*sessionEnd;
        cellID = ones(numPoints(i),1)*i;
        spikes(indices(i)+1:indices(i+1),:) = [tempSpikes cellID];
    end
    cells = spikesToCells(spikes);
end

% groundTruthAssemblies = ad_logicalRand(numCells, sizeAssemblies, numAssemblies);
if randomizeAssemblyMeanSpikes == 2
    lambda = lambda*2;
    lambda(lambda<2) = 2;
end
if exist('assemblyActivationTimestamps','var') && iscell(assemblyActivationTimestamps) && ~isempty(assemblyActivationTimestamps)
    numAActs = zeros(size(assemblyActivationTimestamps));
    for i = 1:length(assemblyActivationTimestamps)
        numAActs(i) = length(assemblyActivationTimestamps{i,1});
    end
else
    numAActs = (abs(ceil(randn(numA,1)*stdAActs + meanAActs)));
end
t = 0;
while any(numAActs < t)
    numAActs(numAActs < t) = randn(sum(numAActs < t),1) * stdAActs2 + (meanAActs);
end

if generateAssemblyIdentity == 1

    groundTruthAssemblies = zeros(numCells, 0);
    numSpikesPerCell = cellfun(@length, cells);
    cellID = (1:numCells)';
    iteration=1; m=0;
    %     listAssemblySize = poissrnd(meanAssemblySize, numA,1);
    %     while any(listAssemblySize < minAssemblySize)
    %         smallerPoints = listAssemblySize<minAssemblySize;
    %         listAssemblySize(smallerPoints) = poissrnd(meanAssemblySize, sum(smallerPoints),1);
    %     end
    listAssemblySize = repelem(sizeAssemblies, numAssemblies);
    listAssemblySize = sort(listAssemblySize,'descend');
    while size(groundTruthAssemblies(:,sum(groundTruthAssemblies)>0),2) < numA
        currAssemblySize= listAssemblySize(iteration);
        if randomizeAssemblyMeanSpikes == 2
            doKeepCell = numSpikesPerCell > (numAActs(iteration)*round(lambda/2)+50);
%             doKeepCell = numSpikesPerCell > (numAActs(iteration)*lambda'+50);
        else
            doKeepCell = numSpikesPerCell > (numAActs(iteration)*assemblyMeanSpikes);
        end
        currSpikesPerCell = numSpikesPerCell(doKeepCell);
        currCellID = cellID(doKeepCell);
        currID = 0;
        while length(unique(currID)) < currAssemblySize
%             currID = randsample(length(currCellID),listAssemblySize(iteration), true, currSpikesPerCell);
            currID = randsample(length(currCellID),listAssemblySize(iteration));
            m = m+1;
            if m== 500
                error(['The dataset cannot take into account the current '...
                    'number and size of chosen assemblies and their '...
                    'activations. Please reduce either of those parameters.'])
            end
        end
        currAssembly = currCellID(currID);
        
        prevNumA = size(groundTruthAssemblies(:,sum(groundTruthAssemblies)>0),2);
        groundTruthAssemblies(currAssembly,iteration) = 1;
        groundTruthAssemblies = unique(groundTruthAssemblies','rows','stable')';
        
        %Remove assemblies contained within others
        v = double(groundTruthAssemblies)' * double(groundTruthAssemblies);
        v(1:length(v)+1:end) = 0;
        z = sum(groundTruthAssemblies);
        p = bsxfun(@eq, v,z);
        groundTruthAssemblies(:,sum(p)>0) = [];
        newNumA = size(groundTruthAssemblies(:,sum(groundTruthAssemblies)>0),2);
        
        if newNumA > prevNumA
            %             numSpikesPerCell(currAssembly) = numSpikesPerCell(currAssembly) - (numAActs(iteration)*lambdaSpikesPerActivation(iteration, currAssembly)');
            if randomizeAssemblyMeanSpikes == 2
                numSpikesPerCell(currAssembly) = numSpikesPerCell(currAssembly) - (numAActs(iteration)*round(lambda(currAssembly)/2));
%                 numSpikesPerCell(currAssembly) = numSpikesPerCell(currAssembly) - (numAActs(iteration)*lambda(currAssembly));
            else
                numSpikesPerCell(currAssembly) = numSpikesPerCell(currAssembly) - (numAActs(iteration)*assemblyMeanSpikes);
            end
            iteration = iteration +1;
            m = 0;
        end
        m = m+1;
        if m== 500
            error(['The dataset cannot take into account the current '...
                'number and size of chosen assemblies and their '...
                'activations. Please reduce either of those parameters.'])
        end
    end
end

groundTruthJitter = zeros(size(groundTruthAssemblies,2),1);
for i = 1:size(groundTruthAssemblies,2)
    groundTruthJitter(i,1) = randn * assembly_stdJitter + assembly_meanJitter;
    while groundTruthJitter(i,1) < assembly_minJitter
        groundTruthJitter(i,1) = randn * assembly_stdJitter + assembly_meanJitter;
    end
end

lambdaSpikesPerActivation = ones(numCells,1)*assemblyMeanSpikes;
cellsUpdated = cells;
newCells = cell(numCells,1);
droppedAssemblies = zeros(size(groundTruthAssemblies,2),1);
assemblyTimestamps = cell(size(groundTruthAssemblies,2),1);
for i = 1:size(groundTruthAssemblies,2)
    currID = find(groundTruthAssemblies(:,i));
    numCurrCells = length(currID);
    
    %% Get timestamps of activations
    
    if ~isempty(assemblyActivationTimestamps)
        if iscell(assemblyActivationTimestamps)
            timestampsAActs = assemblyActivationTimestamps{i,1};
        else
            timestampsAActs = assemblyActivationTimestamps(:,i);
        end

    else
        timestampsAActs = rand(numAActs(i),1)*(sessionEnd-sessionStart) + sessionStart;
    end
    
    assemblyTimestamps{i,1} = timestampsAActs;
    %     jitteredTimestamps = bsxfun(@plus, timestampsAActs, (rand(length(timestampsAActs),numCurrCells)-0.5)*assemblyJitter(i));
    
    %% Find nearest spikes and shift them to new timestamp location +- jitter
    for j = 1:numCurrCells
        if randomizeAssemblyMeanSpikes == 0
            numSpikes = ones(length(timestampsAActs),1) * assemblyMeanSpikes;
        elseif randomizeAssemblyMeanSpikes == 1
            numSpikes = poissrnd(assemblyMeanSpikes, length(timestampsAActs),1);
        elseif randomizeAssemblyMeanSpikes == 2
            numSpikes = ones(length(timestampsAActs),1) * round(lambda(currID(j))/2);
%             numSpikes = poissrnd(lambda(currID(j))/2,length(timestampsAActs),1);
        end
        currTimestamps = repelem(timestampsAActs, numSpikes);
        currTimestamps = currTimestamps + (rand(length(currTimestamps),1)-0.5) * groundTruthJitter(i);
        newCells{currID(j)} = [newCells{currID(j),1}; currTimestamps]; % Add activation timestamps to the new cells
        toRemoverandomizedIDs = randsample(length(cellsUpdated{currID(j)}), length(currTimestamps));
        cellsUpdated{currID(j)}(toRemoverandomizedIDs) = [];
        
%         while ~isempty(toRemoveTimestamps)
%             [~,nearestIDs] = findNearest(cellsUpdated{currID(j)}, toRemoveTimestamps');
%             [uniqueVals,uniqueIndex] = unique(nearestIDs);
%             toRemoveTimestamps(uniqueIndex) = [];
%             cellsUpdated{currID(j)}(uniqueVals) = [];
%         end
    end
    
    if negativeCorrelations == 1
        nonAcells = find(~groundTruthAssemblies(:,i));
        randCell = nonAcells(randi(length(nonAcells),1));
        numSpikes = ones(length(timestampsAActs),1);
        currTimestamps = repelem(timestampsAActs, numSpikes);
        newCells{randCell} = [newCells{randCell,1}; currTimestamps+1]; % Add activation timestamps to the new cells
        toRemoveTimestamps = sort(currTimestamps); % Remove the closest activation to each of the old timestamps
        while ~isempty(toRemoveTimestamps)
            [~,nearestIDs] = findNearest(cellsUpdated{randCell}, toRemoveTimestamps');
            [uniqueVals,uniqueIndex] = unique(nearestIDs);
            toRemoveTimestamps(uniqueIndex) = [];
            cellsUpdated{randCell}(uniqueVals) = [];
        end
    end
end

newSpikes = [cellsToSpikes(newCells); cellsToSpikes(cellsUpdated)];
spikes = sortrows(newSpikes);

%%

if globalActivityVariation == 1 && generateSpikes == 1
    t = spikes(:,1);  % substitute by the timing of your simulated spikes
    muaSpikes = rand(ceil(size(t,1)/200),1)*max(t); % multi-unit activity. Here I've just used rand() of a small number of events (200 times less than 't') within the range of t to have some fluctuations, but it can be substrituted for actual recorded mua if we want
    smooth = 0.01/mode(diff(t)); % 10ms smooth to accentuate the fluctuations a bit (maybe wouldn't be needed with mua of recorded data)

    % Dilate time:
    [h,ht] = hist(muaSpikes,length(muaSpikes)); % histogram of multiunit activity
    h = cumsum(Smooth(h,smooth )); % cumulative multiunit activity
    h = h/max(h)*(max(t)); % scaling it to the same range as 't'. 
    % Now in practice, 'ht' contains uniform ('flat') time, and 'h' contains a transformed version of time that includes the fluctuations. We want to interpolate 't' from flat-time to fluctionations-time:
    [~,u,~] = unique(h); % prevent interpolation error in case of a duplicated value
    spikes(:,1) = interp1(h(u),ht(u),t); % Interpolate time (dilating/shrinking it) so that it reproduces the fluctuations in 'h'.
end


spikes(isnan(spikes(:,1)),:) = [];

