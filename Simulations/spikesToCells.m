function [cells, octNum] = spikesToCells(spikes)
% Transforms spikes, as returned by SetCurrentSession's DATA.spikes, into a
% cell.
%
% Inputs:
% - spikes: in DATA.spikes
%   ---col 1: time.
%   ---col 2: cell number
%   ---col 3: octrode number (optional)
%
% Outputs:
% - cells: cell with the different neurons in columns. Each cell contains
%          the neuron's spike times.
% - octNum: to which octrode each cell belongs
%

numCells = max(spikes(:,2));
cells = cell(numCells,1);
octNum = nan(numCells,1);

if size(spikes,2) == 3
    for i = 1:numCells
        currSpikes = spikes(spikes(:,2) == i, :);
        cells(i,1) = {currSpikes(:,1)};
        octNum(i,1) = mode(currSpikes(:,3));
    end
else
    for i = 1:numCells
        currSpikes = spikes(spikes(:,2) == i, :);
        cells(i,1) = {currSpikes(:,1)};
    end
end
