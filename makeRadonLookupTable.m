function out = makeRadonLookupTable(nSpatialBins, nTemporalBins)
%function out = makeRadonLookupTable(nSpatialBins, nTemporalBins(list))
% Computes a look-up table of the properties of the lines projected by the
% Radon transform. These properties (though not the actual results of the
% transfrom itself) are unique to, and preserved across all matrices with
% of the same size, so its convenient to pre-compute them rather than do so
% for every event.
% Positional inputs:
%       'nSpatialBin' (positive integer scalar): number of spatial bins in
%           the analysis (note that as written, this number is assumed to
%           be constant across all events as is generally the case when
%           analyzing a given replay session/day- if this assumption is
%           broken the code can be run repeatedly with the different
%           spatial-bin-number values)
%       'nTemporalBins': (positive integer vector) containing a list of the
%           unique number of temporal-sub-bins in each event. 
% Output: RadonLookUpTable structure.
%
% Written by Andres Grosmark in 2021


O1 = ones(nSpatialBins, max(nTemporalBins));


theta = deg2rad(0:0.5:179);

allSlopes = 1./tan(theta); 
% this is neccesary due to the odd coordinate system used by the Radon
% transform + the fact that the line we actually care about is perpendicular
% to the one defined by 'theta'

for S = 1:length(nTemporalBins)
    nRadonPoints = 2*ceil(norm([nSpatialBins, nTemporalBins(S)]-floor(([nSpatialBins, nTemporalBins(S)]-1)/2)-1))+3; % this is the equation used by the build in Matlab radon function
    nRadonPoints = nRadonPoints*2; %upscale a bit to get better resolution - this is not strictly neccessary 
    % and is a good candidate to get rid (set 'scale' to 1 instead of '2') if a speed-up is needed
    
    [RO, xp] = radon(O1(:, 1:nTemporalBins(S)), rad2deg(theta), nRadonPoints);
    out{nTemporalBins(S)} = [];
    out{nTemporalBins(S)}.pathLength = RO;
    out{nTemporalBins(S)}.xp = xp;
    out{nTemporalBins(S)}.theta = rad2deg(theta);
    out{nTemporalBins(S)}.nRadonPoints = nRadonPoints;
    
    
    %This is how the radon function defines the center of the image (from
    %which the rest of the polar-coordinate system is laid).
    centerX = floor((nTemporalBins(S) + 1)/2);
    centerY = floor((nSpatialBins + 1)/2);  
    
    out{nTemporalBins(S)}.point1X = NaN(size(RO));
    out{nTemporalBins(S)}.point1Y = NaN(size(RO));
    out{nTemporalBins(S)}.point2X = NaN(size(RO));
    out{nTemporalBins(S)}.point2Y = NaN(size(RO));
    out{nTemporalBins(S)}.size = size(O1(:, 1:nTemporalBins(S)));
    
    out{nTemporalBins(S)}.slope = repmat(allSlopes(:)', [length(xp), 1]);

    for ti = 1:length(theta) %step across the different angles of the lines cast by the radon transform
        for bi = 1:length(xp) %step across the different radial coordinates (relative to the image center)
            [x, y] = pol2cart(theta(ti), xp(bi)); %this gives us one point along the radon line, with this and the slope (calculated above)
            % we can calculate the rest using goold old 8th grade algebra. 
            x2 = centerX + x;
            y2 = centerY - y;
            
            m = allSlopes(ti);
            b = y2 - m*x2; % find the intercept
            
            % calculate the 4 possible points which may intersect the edge
            % of the rectangle defined by the event size:
            topEdge = nSpatialBins + 0.5;
            rightEdge = nTemporalBins(S) + 0.5;
            leftPoint = [0, b];
            rightPoint = [rightEdge, rightEdge*m + b];
            bottomPoint = [-b/m, 0];
            topPoint = [(topEdge - b)/m, topEdge];
            points = [leftPoint; rightPoint; bottomPoint; topPoint];
            % then pick the two points which actually intersect the
            % event-rectangle 
            
            % Coder's note: it feels like there should be a faster way to
            % do this (find the valid points defining the edges of this
            % line) but it did not occur to me, and it doesn't matter much
            % since this look-up table is not calculated on a
            % per-event basis.
            
            k = find(points(:, 1) >= 0 & points(:, 1) <= rightEdge & points(:, 2) >= 0 & points(:, 2) <= topEdge);
            
            if length(k) >= 2
                out{nTemporalBins(S)}.point1X(bi, ti) = points(k(1), 1);
                out{nTemporalBins(S)}.point1Y(bi, ti) = points(k(1), 2);
                out{nTemporalBins(S)}.point2X(bi, ti) = points(k(2), 1);
                out{nTemporalBins(S)}.point2Y(bi, ti) = points(k(2), 2);
            end            
        end
    end
    out{nTemporalBins(S)}.slope(isnan(out{nTemporalBins(S)}.point1X)) = NaN;
    out{nTemporalBins(S)}.pathLengthFromPoints = hypot(out{nTemporalBins(S)}.point2X - out{nTemporalBins(S)}.point1X, ...
        out{nTemporalBins(S)}.point2Y - out{nTemporalBins(S)}.point1Y);
    out{nTemporalBins(S)}.spaceOffset = out{nTemporalBins(S)}.point2Y - out{nTemporalBins(S)}.point1Y;
    out{nTemporalBins(S)}.tempOffset = out{nTemporalBins(S)}.point2X - out{nTemporalBins(S)}.point1X;
    out{nTemporalBins(S)}.spaceOffsetRound = round(out{nTemporalBins(S)}.spaceOffset);
    out{nTemporalBins(S)}.tempOffsetRound = round(out{nTemporalBins(S)}.tempOffset);
    out{nTemporalBins(S)}.tempOffsetRoundPerc = 100*abs(floor(out{nTemporalBins(S)}.tempOffset))/nTemporalBins(S);
end