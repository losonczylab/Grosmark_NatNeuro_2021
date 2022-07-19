function synthEvents = createSyntheticEvents(varargin)
% This function generates 1) synthetic place fields, and 2) synthetic population 
% activity events for demonstrating replay analysis. By default the current
% function generates 5 events characterizing distinct types of replay or
% noise signals: 1) Linear Forward Replay, 2) Reverse Linear Replay, 3)
% Forward Circular (i.e. edge-spanning) Replay, 4) Reverse Circular Replay,
% 5) Unstructured (noise) Activity
% Input the following optional paramaters may be inputed as string-value pairs:
%   'nSpatialBins' (default = 100), number of total spatial bins for generating
%           synthetic place fields 
%   'nPlaceCells' (default = 30) number of synthetic place cells. 
%   'peakRate' (deault = 10) peak place field firing rate of cells
%   'nTimeBinsPerEvent' (default = 15) number of (non-overlapping) bins 
%           within each synthetic population event
%   'eventBinDuration' (default 0.03) duration (in seconds) of the bins
%           within each event. This is needed because rate (rather than count
%           number) is used for Bayesian decoding. 
%
% Outputs: a structure with the following fields:
%       'params': a list of the paramaters used for generating the
%           synthetic events and synthetic place fields.
%       'rateByPositionMatrix': an [nCell X nSpatialBin] matrix of firing
%           rates (i.e. synthetic place fields)
%       'eventFiringRateMatrix': [nTimeBinsPerEvent*nEvents(def=5), nCell]
%           matrix of within-event by-time-bin firing rates
%       'eventFiringRateMatrixID' = [nTimeBinsPerEvent*nEvents(def=5), 1] 
%           column vector of the unique numeric ids indicating which unique
%           event each of the bins in 'eventFiringRateMatrix' correspond
%           to.
%       'eventType': a cell array of labels indicating which type of replay
%           phenomenon is being simulated by each synthetic events 
%       'eventTypeMultiLine': same as 'eventType' but split up into
%           2 lines to make them easier to use as panel titles. 
%
% Written by Andres Grosmark in 2021

nSpatialBins = 100;
placeFieldGaussWidth = 0.05;
nPlaceCells = 50;
peakRate = 10;
nTimeBinsPerEvent = 15;
eventBinDuration = 0.03;
for i = 1:(length(varargin) - 1)
    if ischar(varargin{i})
        switch lower(varargin{i})
            case lower('nSpatialBins')
                nSpatialBins = varargin{i + 1};
            case lower('placeFieldGaussWidth')
                placeFieldGaussWidth = varargin{i};
            case lower('nPlaceCells')
                placeFieldGaussWidth = nPlaceCells;
            case lower('peakRate')
                peakRate = varargin{i + 1};
            case lower('nTimeBinsPerEvent')
                nTimeBinsPerEvent = varargin{i + 1};
            case lower('eventBinDuration')
                eventBinDuration = varargin{i + 1};
        end
    end
end

synthEvents = [];
synthEvents.params.nSpatialBins = nSpatialBins;
synthEvents.params.placeFieldGaussWidth = placeFieldGaussWidth;
synthEvents.params.nPlaceCells = nPlaceCells;
synthEvents.params.peakRate = peakRate;
synthEvents.params.nTimeBinsPerEvent = nTimeBinsPerEvent;
synthEvents.params.eventBinDuration = eventBinDuration;



%% Create synthetic place fields
if mod(nSpatialBins, 2) == 0 % this is to make sure our gaussian have only 1 peak value
    placeFieldGauss = fspecial('Gaussian', [nSpatialBins + 1, 1], placeFieldGaussWidth*nSpatialBins);
    placeFieldGauss = placeFieldGauss(1:(end - 1));
else
    placeFieldGauss = fspecial('Gaussian', [nSpatialBins, 1], placeFieldGaussWidth*nSpatialBins);
end

% normalize our gaussian to the peak rate
placeFieldGauss = peakRate*(placeFieldGauss/max(placeFieldGauss));
[~, peakLoc] = max(placeFieldGauss);
placeFieldGauss = circshift(placeFieldGauss, -(peakLoc - 1));


pfPeakLocations = round(linspace(0, nSpatialBins - 1, nPlaceCells));

synthEvents.rateByPositionMatrix = [];
for i = 1:nPlaceCells
    synthEvents.rateByPositionMatrix = [synthEvents.rateByPositionMatrix; circshift(placeFieldGauss, pfPeakLocations(i))'];
end


%% Create 5 types of synthetic population events

synthEvents.eventFiringRateMatrix = [];
synthEvents.eventFiringRateMatrixID = [];

synthEvents.eventType = {};
synthEvents.eventTypeMultiLine = {};
cellBinInEvent = ceil(((pfPeakLocations + 1)/nSpatialBins)*nTimeBinsPerEvent);

for candidateEvent = 1:5
    eventFR = zeros(nTimeBinsPerEvent, nPlaceCells);
    switch candidateEvent
        case 1
            % simple forward linear replay
            for i = 1:length(cellBinInEvent)
                eventFR(cellBinInEvent(i), i) = 1;
            end
            eventFR = eventFR/eventBinDuration;
            synthEvents.eventFiringRateMatrix = [synthEvents.eventFiringRateMatrix; eventFR];
            synthEvents.eventFiringRateMatrixID = [synthEvents.eventFiringRateMatrixID; repmat(candidateEvent, [size(eventFR, 1), 1])];
            synthEvents.eventType{end + 1}  = 'Linear Forward Replay';
            synthEvents.eventTypeMultiLine{end + 1} = {'Linear', 'Forward Replay'};
            
        case 2
            % Linear Reverse Replay -- same as case 1 above but with the
            % time axis (dim = 1) flippled. 
            for i = 1:length(cellBinInEvent)
                eventFR(cellBinInEvent(i), i) = 1;
            end
            eventFR = eventFR/eventBinDuration;
            eventFR = flipud(eventFR); %flip forward event to make reverse event
            synthEvents.eventFiringRateMatrix = [synthEvents.eventFiringRateMatrix; eventFR];
            synthEvents.eventFiringRateMatrixID = [synthEvents.eventFiringRateMatrixID; repmat(candidateEvent, [size(eventFR, 1), 1])];
            synthEvents.eventType{end + 1}  = 'Linear Reverse Replay';
            synthEvents.eventTypeMultiLine{end + 1} = {'Linear', 'Reverse Replay'};
            
        case 3
            % Circular Forward Replay -- here we circularly offset the linear event above such that
            % it contains circular but not linear replay content
            cellBinInEventCircShift = circshift(cellBinInEvent(:),  round(nPlaceCells/4));
            for i = 1:length(cellBinInEventCircShift)
                eventFR(cellBinInEventCircShift(i), i) = 1;
            end
            eventFR = eventFR/eventBinDuration;
            synthEvents.eventFiringRateMatrix = [synthEvents.eventFiringRateMatrix; eventFR];
            synthEvents.eventFiringRateMatrixID = [synthEvents.eventFiringRateMatrixID; repmat(candidateEvent, [size(eventFR, 1), 1])];
            synthEvents.eventType{end + 1}  = 'Circular Forward Replay';
            synthEvents.eventTypeMultiLine{end + 1} = {'Circular', 'Forward Replay'};
        case 4
            % Circular Reverse Replay -- same as case 3 above but with the
            % time axis (dim = 1) flippled. 
            cellBinInEventCircShift = circshift(cellBinInEvent(:),  round(nPlaceCells/4));
            for i = 1:length(cellBinInEventCircShift)
                eventFR(cellBinInEventCircShift(i), i) = 1;
            end
            eventFR = eventFR/eventBinDuration;
            eventFR = flipud(eventFR); %flip forward event to make reverse event
            synthEvents.eventFiringRateMatrix = [synthEvents.eventFiringRateMatrix; eventFR];
            synthEvents.eventFiringRateMatrixID = [synthEvents.eventFiringRateMatrixID; repmat(candidateEvent, [size(eventFR, 1), 1])];
            synthEvents.eventType{end + 1}  = 'Circular Reverse Replay';
            synthEvents.eventTypeMultiLine{end + 1} = {'Circular', 'Reverse Replay'};
            
        case 5
            % Random acitvity - here we randomly reassign each group of
            % cells coding for the same time-bin to a new time bin. This
            % method preserves the instantaneous correlation structure and
            % therefeore preserves the 'peakiness' of the resulting
            % Bayesian reconstruction. Note that this is more straight-forward
            % to do by permuting the Bayesian posteriors themselves (see
            % 'time-swap' shuffle below which has the same result) but I
            % chose to do it this way to keep these examples consistent.
            cellBinInEventShuff = cellBinInEvent;
            uBins = unique(cellBinInEvent);
            uBinsShuff = uBins(randperm(length(uBins)));
            for i = 1:length(uBins)
                cellBinInEventShuff(cellBinInEvent == uBins(i)) = uBinsShuff(i);
            end
            
            for i = 1:length(cellBinInEventShuff)
                eventFR(cellBinInEventShuff(i), i) = 1;
            end
            eventFR = eventFR/eventBinDuration;
            synthEvents.eventFiringRateMatrix = [synthEvents.eventFiringRateMatrix; eventFR];
            synthEvents.eventFiringRateMatrixID = [synthEvents.eventFiringRateMatrixID; repmat(candidateEvent, [size(eventFR, 1), 1])];
            synthEvents.eventType{end + 1} = 'Unstructured Activity';
            synthEvents.eventTypeMultiLine{end + 1} = {'Unstructured', 'Activity'};
    end
end