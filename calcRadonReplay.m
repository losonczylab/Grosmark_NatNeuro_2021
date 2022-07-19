function [posMean, id, slope, pathLength, pointsXYXY, MaxID] = calcRadonReplay(posteriors, binIDs, varargin)
% [posMean, id, slope, pathLength, pointsXYXY] = calcRadonReplay(posteriors, binIDs, radonLookUpTable(optional))
% Calculates the 'line-casting' (Radon transform) replay properties across
% candidate events. Briefly, this equates to densly sampling lines at
% with different slopes and intercepts to see which one of these lines
% intersects the most average posterior probability. In practice it is
% neccesary to restrict the search to those crossing at least a majority of
% the temporal bins otherwise lines crossing a single corner-bin tend to be
% selected (see the 'minNBinPerc' hard-coded variable below which is set to
% 100% by default). 
%
% Positional inputs:
%   'posteriors': [nTemporalBin X nSpatialBin] matrix of Posterior
%           probability estimates
%   'binIDs': ['nSpatialBin X 1] column vector of unique (numerical id's)
%           indicating which event each of the temporal bins (rows) in
%           'posteriors' belongs to. 
%   'radonLookupTable' (optional): pre-calculated look-up table with the
%            properties of the densely sampled lines. If not provided,
%            'makeRadonLookupTable' is called within the function instead.
%
%   Outputs:
%       'posMean': [nEvent X 1]average posterior probability intersected by the
%               best-fit line
%       'id': [nEvent X 1] the unique numerical identifier of each event (see
%               'binIDs')
%       'slope': [nEvent X 1] the slope (in bin-units) of the best-fit line
%       'pathLength': [nEvent X 1] the total length (within the event) of
%              the best fit line
%       'pointsXYXY': [nEvent X 4] an [x(1), y(1), x(2), y(2)] set of
%              starting and ending points within the event (in bin-units)
%              for the best fit line
%       'MaxID': the id of the best-fit line within the Radon-lookup table
%               not currently used outside this function but could be
%               useful for debugging. 
%
%  See: Davidson TJ, Kloosterman F, Wilson MA. Hippocampal replay of
%       extended experience. Neuron. 2009 Aug 27;63(4):497-507. For
%       original implementation.
%
% Written by Andres Grosmark in 2021


minSpatialDisp = 0; % minimum spatial displacement of valid lines (in bins)
minNBinPerc = 100; %minimum percentage of temporal bins that valid lines 
% must cross. 

id = unique(binIDs);
hID = histc(binIDs, id);
uNS = unique(hID);

if ~isempty(varargin)
    radonLookupTable = varargin{1};
else
    radonLookupTable = makeRadonLookupTable(size(posteriors, 1), uNS);
end

posA = [];


posMean = NaN(length(id), 1);
slope = posMean;
pathLength = slope;
pointsXYXY = NaN(length(id), 4);
MaxID = slope;
for U = 1:length(uNS)
    uh = find(hID == uNS(U));
    minPLength = sqrt((floor(uNS(U).*(minNBinPerc/100)) + 1)^2 + minSpatialDisp^2); 
    %% find all the valid lines:
    goodLines = abs(radonLookupTable{uNS(U)}.spaceOffset) >= minSpatialDisp & ...
        abs(radonLookupTable{uNS(U)}.tempOffsetRoundPerc) >= minNBinPerc & ...
        radonLookupTable{uNS(U)}.pathLength >= minPLength;
    theta = radonLookupTable{uNS(U)}.theta;
    nRadonPoints = radonLookupTable{uNS(U)}.nRadonPoints;
    for i = 1:length(uh)
        
        radTrans = radon(posteriors(:, binIDs == id(uh(i))), theta, nRadonPoints);
        
        radTrans(~goodLines) = 0;
        
        % The Radon transform returns a sum, we can calculate a mean by
        % the dividing by the path-length of that line:
        radTrans = radTrans./radonLookupTable{uNS(U)}.pathLength;
        
        [posMean(uh(i)), L] = max(radTrans(:));
        slope(uh(i)) = radonLookupTable{uNS(U)}.slope(L);
        MaxID(uh(i)) = L; 
        pathLength(uh(i)) = radonLookupTable{uNS(U)}.pathLength(L);
        pointsXYXY(uh(i), :) = [radonLookupTable{uNS(U)}.point1X(L), radonLookupTable{uNS(U)}.point1Y(L), ...
            radonLookupTable{uNS(U)}.point2X(L), radonLookupTable{uNS(U)}.point2Y(L)];
    end
end

k = find(pointsXYXY(:, 3) < pointsXYXY(:, 1));
pointsXYXY(k, :) = [pointsXYXY(k, 3:4), pointsXYXY(k, 1:2)];