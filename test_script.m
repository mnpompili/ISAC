
timescale = 0.1; % e.g. 100ms assemblies
spikes; % expecting a 2-column matrix [timestamps, IDs]
intervals = [0 Inf]; % optionally, intervals of interest used in detecting the assemblies (e.g., running periods, ripples, etc.)
% Periods outside of these intervals will be ignored in assembly detection
% (but we will detect the activation of the detected assemblies over all
% times, including outside of the "intervals" of interest).

% Detect ISAC assemblies
d = [0; cumsum(diff(intervals,1,2))];
restricted = cell2mat(arrayfun(@(i) ...
    [spikes(spikes(:,1)>=intervals(i,1) & spikes(:,1)<=intervals(i,2),1) ...
     spikes(spikes(:,1)>=intervals(i,1) & spikes(:,1)<=intervals(i,2),2)], ...
    (1:size(intervals,1))','UniformOutput',0));

assemblies = ISAC(restricted,timescale);

% Find the activations of these ISAC assemblies
activations = FindActivations(spikes,assemblies,timescale);


















