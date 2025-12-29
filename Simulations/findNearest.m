function [nearestValues, nearestIdx] = findNearest(mainMatrix, matchVector, mainColumn)
% Gets you the nearest values to matchVector in mainMatrix (and its
% corresponding row). 
%
% Highly optimized function, uses histogram edges to determine proximity.
%
% If mainMatrix was a matrix, allows you to specify a specific column with
% which to match matchVector. All columns corresponding to nearest value
% are returned in this case.
%
% Inputs:
% - mainMatrix: Vector or matrix to which elements will be matched in a
% certain column.
% - matchVector: Values to be matched.
% - mainColumn: Column of mainMatrix corresponding to matchVector.
%
% Output:
% - nearestValues: values of mainMatrix's mainColumn nearest to
% matchVector.
%


%% Initialize variables
if ~exist('mainColumn', 'var')
    mainColumn = 1;
end

colMatch = mainMatrix(:,mainColumn);

%% Find nearest points using a histogram-based method
edges = [-Inf, mean([colMatch(2:end) colMatch(1:end-1)],2)', +Inf];
nearestIdx = discretize(matchVector, edges);
nearestValues = mainMatrix(nearestIdx,:);
