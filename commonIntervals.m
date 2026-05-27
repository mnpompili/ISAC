function intervals = commonIntervals(spikes, windowSize, numCommon, numCells, doOverlap, type, reverse)

% Get intervals where more than numCommon cells overlap.
% 
% outputIntervals = commonIntervals(spikes, windowSize, numCommon, numCells)
% 
% INPUTs: 
% - spikes: N x 2 matrix. The first column is times, second column
%   is spike identity. 
% - windowSize: size of the window in which to retain overlapping elements
% - numCommon: number of overlapping elements to retain an interval.
% - numCells: number of neurons.
% - doOverlap: set to 1 = the last spike of activation A is the first
% spike of activation B.
% - type:
% ---'min': searches for the smallest interval when numCommon cells are
% active
% ---'max': searches for the largest interval (smaller than windowSize)
% where numCommon cells are active.
% - reverse: if 1: searches for activations starting from the end of the
% session
%
% OUTPUT:
% - outputIntervals: M x 2 matrix of the different intervals which have
%   numCommon overlapping elements.
%
% Copyright (C) 2020-2026 by Gabriel Makdah
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.


if ~exist('doOverlap','var') || isempty(doOverlap)
    doOverlap = 0;
end

if ~exist('type','var') || isempty(type)
    type = 'min';
end

if ~exist('reverse','var') || isempty(reverse)
    reverse = 0;
end

if reverse == 1
    spikes(:,1) = -spikes(:,1);
    spikes = flipud(spikes);
    savedMin = spikes(1,1);
    spikes(:,1) = spikes(:,1) - savedMin;
end

switch type
    case 'min'
        intervals = commonIntervals_min(spikes, windowSize, numCommon, numCells, doOverlap);
    case 'max'
        intervals = commonIntervals_max(spikes, windowSize, numCommon, numCells, doOverlap);
end

if reverse == 1
    intervals = -rot90(savedMin+intervals,2);
end
