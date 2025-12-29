function spikes = cellsToSpikes(cells, k)
% Transform a cell to an M x 2 spikes matrix with the first column spikes
% times, and second column the identity of the spikes.
%
% INPUTs:
% - cells: N x 1 cell. Each cell contains the spikes times of a neuron. The
%   cell's index will subsequently be used for the identity of the spikes
%   (second column of the output).
% - k: optinal scalar, or 1 x m vector of numbers to add to the spikes matrix
%
% OUPUT:
% - spikes: Matrix where the first column is spike times, and second column
% the spike identity. Subsequent column represent K.
%   

%%
if ~exist('k','var') || isempty(k)
    idxCell = cell(size(cells));
    for i = 1:length(cells)
        idxCell{i,1} = ones(length(cells{i,1}),1)*i;
    end

    spikes = [cell2mat(cells) cell2mat(idxCell)];
    spikes = sortrows(spikes);
else
    idxCell = cell(size(cells));
    for i = 1:length(cells)
        idxCell{i,1} = ones(size(cells{i,1},1),1)*[i k];
    end
    
    spikes = [cell2mat(cells) cell2mat(idxCell)];
    spikes = sortrows(spikes);
end