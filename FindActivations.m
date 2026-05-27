function activations = FindActivations(spikes,assemblies,windowSize)

% FindActivations finds the activation events of detected neuron assemblies.
%
%   activations = FindActivations(spikes,assemblies,windowSize)
%
%   This function identifies the time intervals during which each detected
%   assembly was active. An assembly activation is defined as an event in
%   which all neurons belonging to a given assembly fired within the
%   specified temporal window.
%
%   REQUIRED INPUTS
%   spikes         - a two-column [timestamp, unitID] matrix containing
%                    the list of spikes for each unit
%   assemblies     - a binary matrix where each row corresponds to an
%                    assembly and each column corresponds to a neuron.
%                    A value of 1 indicates that the neuron belongs to the
%                    assembly.
%   windowSize     - the maximum allowed temporal gap between spikes for
%                    them to be considered part of the same assembly
%                    activation event.
%
%
%   OUTPUT
%   activations    - a cell array containing the activation events for
%                    each assembly. Each cell contains a matrix where rows
%                    correspond to activation events and columns contain:
%
%                    [start stop] timestamps of the detected activation
%
% Copyright (C) 2026 by Ralitsa Todorova & Gabriel Makdah
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


% Find the activations of each assemblies
for i=1:size(assemblies,1)
    ok = ismember(spikes(:,2),find(assemblies(i,:)));
    code = cumsum(assemblies(i,:))';
    activations{i,1} = commonIntervals([spikes(ok,1) code(spikes(ok,2))], windowSize, sum(assemblies(i,:)), sum(assemblies(i,:)));
end

