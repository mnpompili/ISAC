% One example data set (spikes, groundTruthAssemblies) per simulation of the
% simulation table, in the order of the table, made with
% script_createPoissonSimulations.m and script_createShuffledSimulations.m
% (rngCounter as in runSimulations_base.m).
% Before running:
%  - in both scripts, comment out "whichSimulation = 1;" and "rngCounter = 1;"
%    (otherwise every simulation is simulation 1)
%  - in script_createPoissonSimulations.m, comment out the old melody branch
%    (lines 126-131, "elseif whichSimulation >= (length(numAvar)+..."), which
%    stops simulations 33-39 with an error
% Run it from this folder, so that the corrected melody_getTimestamps.m next
% to it is used (MATLAB uses the current folder first).

numCycle = 1;
listOfSims = [1 2 47 3 48 4 49 5, ...                          % table 1-8:   number of assemblies
              6 33 7 42 8 43 9 44 45, ...                      % table 9-17:  assembly size
              10:19, ...                                       % table 18-27: session length
              20 21 22 23 34 24 35 25 36 37 38 26 39 27, ...   % table 28-41: number of activations
              62:-1:51];                                       % table 42-53: melody
disp(['melody_getTimestamps used: ' which('melody_getTimestamps')])
disp('Starting')
for background = {'poisson', 'peyrache'}
    simulations = struct('simulation', {}, 'whichSimulation', {}, 'rngCounter', {}, 'spikes', {}, 'groundTruthAssemblies', {});
    for tableNumber = 1:length(listOfSims)
        whichSimulation = listOfSims(tableNumber);
        rngCounter = whichSimulation*1000 + numCycle;
        simulations(tableNumber).simulation = tableNumber;
        simulations(tableNumber).whichSimulation = whichSimulation;
        simulations(tableNumber).rngCounter = rngCounter; % saved here: the scripts change rngCounter
        if strcmp(background{1}, 'poisson')
            script_createPoissonSimulations
        else
            script_createShuffledSimulations
        end
        if whichSimulation ~= listOfSims(tableNumber) || rngCounter < simulations(tableNumber).rngCounter
            error('The script reset whichSimulation or rngCounter: comment out "whichSimulation = 1;" and "rngCounter = 1;" in it')
        end
        if retry ~= 0, error('Simulation %d did not complete', whichSimulation); end
        simulations(tableNumber).spikes = spikes;
        simulations(tableNumber).groundTruthAssemblies = groundTruthAssemblies;
        fprintf('%s | table %2d | whichSimulation %2d | %7d spikes | %d cells | %3d assemblies\n', background{1}, ...
            tableNumber, whichSimulation, size(spikes,1), size(groundTruthAssemblies,1), size(groundTruthAssemblies,2));
    end
    if strcmp(background{1}, 'poisson'), poisson = simulations; else, peyrache = simulations; end
end
save('ISAC_example_simulations.mat', 'poisson', 'peyrache', '-v7');
disp('Saved ISAC_example_simulations.mat')
