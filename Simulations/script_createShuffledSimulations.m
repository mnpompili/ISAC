
% Each simulations will produce simulated assembly activations in the
% background activity of neurons with realistic statistics produced by
% shuffling prefrontal recordings

% The following simulation parameters are considered default parameters:
% session = 3600 s (1 hour)
% numCells = 50 % number of simulated cells
% numAssemblies = 20 % number of simulated assemblies
% sizeAssemblies = 5 % number of assembly members per each assembly
% meanAActs = 120 % mean number of activations for each assembly
%
% Particular simulation parameters are modified by changing the value of the
% "whichSimulation" variable according to the following table of effects.
% The table of effects indicates, for each value of the "whichSimulation"
% variable, which parameter will be modified from the default (other
% parameters stay the default value).
%
% ----whichSimulation = [1:5, 47:50]: varying the number of simulated assemblies----
% 1 | numAssemblies = 10
% 47 | numAssemblies = 35
% 3 | numAssemblies = 50 
% 48 | numAssemblies = 75
% 4 | numAssemblies = 100 
% 49 | numAssemblies = 150
% 5 | numAssemblies = 200
% 50 | numAssemblies = 300
%
% ----whichSimulation = [6:9, 33, 42:46]: varying simulated assembly size while keeping the number of member neurons fixed----
% 6 | sizeAssemblies = 3; numAssemblies = 67
% 33| sizeAssemblies = 4; numAssemblies = 33
% 7 | sizeAssemblies = 5; numAssemblies = 20
% 8 | sizeAssemblies = 7; numAssemblies = 10
% 9 | sizeAssemblies = 10; numAssemblies = 4 
% 42 | sizeAssemblies = 6; numAssemblies = 13
% 43 | sizeAssemblies = 8; numAssemblies = 7
% 44 | sizeAssemblies = 12; numAssemblies = 3
% 45 | sizeAssemblies = 15; numAssemblies = 2
% 46 | sizeAssemblies = 20; numAssemblies = 1
%
% ----whichSimulation = [10:19]: varying simulated recording duration while keeping the assembly activation rate fixed----
% 10 | sessionEnd = 10 minutes (600 s); meanAActs = 10 
% 11 | sessionEnd = 20 minutes (1200 s); meanAActs = 20 
% 12 | sessionEnd = 30 minutes (1800 s); meanAActs = 30 
% 13 | sessionEnd = 40 minutes (2400 s); meanAActs = 40 
% 14 | sessionEnd = 50 minutes (3000 s); meanAActs = 50 
% 15 | sessionEnd = meanAActs = 60 
% 16 | sessionEnd = 70 minutes (4200 s); meanAActs = 70 
% 17 | sessionEnd = 80 minutes (4800 s); meanAActs = 80 
% 18 | sessionEnd = 90 minutes (5400 s); meanAActs = 90 
% 19 | sessionEnd = 100 minutes (6000 s); meanAActs = 100 
% 
% ----whichSimulation = [20:27, 35:39]: varying the mean number of assembly activations (as a fraction of recording duration)----
% 20 | 0.05 activations/min i.e. meanAActs = 3
% 21 | 0.1 activations/min i.e. meanAActs = 6
% 22 | 0.2 activations/min i.e. meanAActs = 12
% 23 | 0.5 activations/min i.e. meanAActs = 30
% 34 | 0.75 activations/min i.e. meanAActs = 45
% 24 | 1 activations/min i.e. meanAActs = 60
% 35 | 1.5 activations/min i.e. meanAActs = 90
% 25 | 2 activations/min i.e. meanAActs = 120
% 36 | 2.5 activations/min i.e. meanAActs = 150
% 37 | 3 activations/min i.e. meanAActs = 180
% 38 | 4 activations/min i.e. meanAActs = 240
% 26 | 5 activations/min i.e. meanAActs = 300
% 39 | 7 activations/min i.e. meanAActs = 420
% 27 | 10 activations/min i.e. meanAActs = 600
%
% ----whichSimulation = 40-41: short recordings----
% 40 | sessionEnd = 2 minutes (120 s); meanAActs = round(120/60) = 2
% 41 | sessionEnd = 5 minutes (300 s); meanAActs = round(300/60) = 5
%
% ----whichSimulation >= 51: a particular simulation in which assembly timestamps driven by melody----
% consult the melody_getTimestamps script


whichSimulation = 1; % select which simulation to run (consult table above)
rngCounter = 1; % set the random number generator for easy replicability

rng(rngCounter)

filepath = fileparts(mfilename('fullpath'));

data = load(fullfile(filepath,'Peyrache2009_spikes.mat'),'spikes'); spikes = data.spikes;
cells = spikesToCells(spikes);
spkRate = cellfun(@length,cells)/(spikes(end,1));
cells(spkRate == 0) = [];

spikes = cellsToSpikes(cells);
for i = 1:2
    spikes = [spikes; spikes(:,1)+max(spikes(:,1)) spikes(:,2)];
end
spikes = ad_shuffleID(spikes);
cells = spikesToCells(spikes);
% cells = ad_shuffleISI(spikesToCells(spikes));
[~,cellID] = sort(cellfun(@length,cells),'descend');


sessionStart = 0;
sessionEnd = 3600;
numCells = 50;
sizeAssemblies = 5;
% disp('Warning, running with 10 assemblies instead of 20')
numAssemblies =  20;
assemblyMeanSpikes = 2;

meanAActs = 120; % Average number of times assemblies are active
stdAActs = 0; % SD of the first run to compute the number of times assemblies are active
stdAActs2 = 0; % SD of second run if cells have negative firing rates
generateAssemblyIdentity = 1; % 1: generate assemblies || 0: if you supply groundTruthAssemblies
randomizeAssemblyMeanSpikes = 2; %  0: exactly assemblyMeanSpikes every assembly activation || 1: poisson point process with lambda = assemblyMeanSpikes || 2: poisson point process with lambda = cellFiringRate/ 2
globalActivityVariation = 0; % Set to 1 to include global variations in cell firing rates
setMeanCellActivity = 0; % set to 0 if you're supplying cell firing rates manually (in cells)
generateSpikes = 0; % set to 0 if you're supplying spikes manually
negativeCorrelations = 0; % set to 1 to include random negative correlations within the dataset
assemblyActivationTimestamps = {};

assembly_meanJitter = 0.03;
assembly_stdJitter = 0;
assembly_minJitter = 0;

firingRate_mean = 2;
firingRate_jitter = 0;
firingRate_jitter2 = 0;

numAvar = [10 20 50 100 200];
sizeAvar = [3 5 7 10];
lengthSessVar = 10:10:100;
numActVar = [0.05 0.1 0.2 0.5 1 2 5 10];
numMelodyActs = [1:12 0.75 0.5];
totalCombs = numAssemblies*(factorial(sizeAssemblies)/factorial(sizeAssemblies-2)*factorial(2));

if whichSimulation >= 1 && whichSimulation < (length(numAvar)+1)
        numAssemblies = numAvar(whichSimulation);
elseif whichSimulation >= (length(numAvar)+1) &&...
        whichSimulation < (length(numAvar)+length(sizeAvar)+1)
    sizeAssemblies = sizeAvar(whichSimulation-length(numAvar));
    numAssemblies = round(totalCombs/(factorial(sizeAssemblies)/factorial(sizeAssemblies-2)*factorial(2)));
elseif whichSimulation >= (length(numAvar)+length(sizeAvar)+1) &&...
        whichSimulation < (length(numAvar)+length(sizeAvar)+length(lengthSessVar)+1)
    sessionEnd = lengthSessVar(whichSimulation-length(numAvar)-length(sizeAvar)) * 60;
    meanAActs = round(sessionEnd/60);
elseif whichSimulation >= (length(numAvar)+length(sizeAvar)+length(lengthSessVar)+1) &&...
        whichSimulation < (length(numAvar)+length(sizeAvar)+length(lengthSessVar)+length(numActVar)+1)
    meanAActs = round(numActVar(whichSimulation-length(numAvar)-length(sizeAvar)-length(lengthSessVar))*sessionEnd/60);
elseif whichSimulation == 33
    sizeAssemblies = 4;
    numAssemblies = round(totalCombs/(factorial(sizeAssemblies)/factorial(sizeAssemblies-2)*factorial(2)));
elseif whichSimulation == 34
    meanAActs = round(0.75*sessionEnd/60);
    elseif whichSimulation == 35
    meanAActs = round(1.5*sessionEnd/60);
    elseif whichSimulation == 36
    meanAActs = round(2.5*sessionEnd/60);
    elseif whichSimulation == 37
    meanAActs = round(3*sessionEnd/60);
    elseif whichSimulation == 38
    meanAActs = round(4*sessionEnd/60);
    elseif whichSimulation == 39
    meanAActs = round(7*sessionEnd/60);
elseif whichSimulation == 40
    sessionEnd = 2 * 60;
    meanAActs = round(sessionEnd/60);
    elseif whichSimulation == 41
    sessionEnd = 5 * 60;
    meanAActs = round(sessionEnd/60);
elseif whichSimulation == 42
    sizeAssemblies = 6;
    numAssemblies = round(totalCombs/(factorial(sizeAssemblies)/factorial(sizeAssemblies-2)*factorial(2)));
    elseif whichSimulation == 43
    sizeAssemblies = 8;
    numAssemblies = round(totalCombs/(factorial(sizeAssemblies)/factorial(sizeAssemblies-2)*factorial(2)));
    elseif whichSimulation == 44
    sizeAssemblies = 12;
    numAssemblies = round(totalCombs/(factorial(sizeAssemblies)/factorial(sizeAssemblies-2)*factorial(2)));
    elseif whichSimulation == 45
    sizeAssemblies = 15;
    numAssemblies = round(totalCombs/(factorial(sizeAssemblies)/factorial(sizeAssemblies-2)*factorial(2)));
    elseif whichSimulation == 46
    sizeAssemblies = 20;
    numAssemblies = round(totalCombs/(factorial(sizeAssemblies)/factorial(sizeAssemblies-2)*factorial(2)));
elseif whichSimulation == 47
    numAssemblies = 35;
    elseif whichSimulation == 48
    numAssemblies = 75;
    elseif whichSimulation == 49
    numAssemblies = 150;
    elseif whichSimulation == 50
    numAssemblies = 300;
elseif whichSimulation >= 51
    assemblyActivationTimestamps = melody_getTimestamps(numMelodyActs(whichSimulation-50));
    numAssemblies = size(assemblyActivationTimestamps,1);
    % sessionEnd = 600;
    sessionEnd = 1200;
    numCells = 20;
else
    error('Value exceeds bound')
end

cells = cells(cellID(1:numCells));
spikes = cellsToSpikes(cells);
retry = 1;
tooMuch = 0;
while retry == 1 && tooMuch < 20
    try
        tooMuch = tooMuch+1;
        rngCounter = rngCounter +100000;
        rng(rngCounter)
        simgetDet
        retry = 0;
    end
    if tooMuch == 100
        error('Could not complete the simulation')
    end
end
