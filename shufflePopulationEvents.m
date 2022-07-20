function [shuffledEvents, newRateMatrix] = shufflePopulationEvents(shuffleType, varargin)
% function [shuffledEvents, newRateMatrix] = shufflePopulationEvents(shuffleType, varargin)
% Creates shuffled (or synthetic versions) of population activity events. There are
% four types of shuffle implemented (the choice must be inputed as a string as the first
% positional input argument:
%       1) 'timeBinPermutation': permutes (resamples without replacement)
%                   Bayesian posterior estimates across the time bins within an event
%                    -Required string-value pair inputs: 'bayesPostProb' and 'binIDs'
%
%       2) 'withinEventActiveCellIDPermutation': shuffles the IDs of only
%                   those subset of cells active within an event, and then re-performs 
%                   the Bayesian decoding on these shuffled ids
%                   -Required string-value pair inputs:
%                   'eventFiringRateMatrix', 'binIDs',
%                   'rateByPositionMatrix', 'eventBinDuration'
%
%       3) 'crossEventPoissonSurrogate': creates a new set of synthetic
%                   randomly generated using Poisson distributions whos
%                   rate parameters match those of the observed cells
%                   across events. This is performed in an event (rather
%                   than rate) based way by first multiplying the rate by the
%                   'eventBinDuration' scalar. Note: I do not think this is
%                   a particularly good shuffle because spiking is not
%                   strictly speaking Poisson in general, and particularly
%                   not so within events which are typically artifically
%                   selected to have high firing rates (see Grosmark &
%                   Buzsaki 2016, Fig S8). However, since this was the
%                   method used in an influential study (Silva et al. 2015)
%                   the method is included here.
%                   -Required string-value pair inputs:
%                   'eventFiringRateMatrix', 'binIDs',
%                   'rateByPositionMatrix', 'eventBinDuration'
%
%       4) 'circularPlaceFieldShuffle': randomly circularly rotates the
%                   place fields and recomputes the Bayesian posteror
%                   decoding (see Grosmark & Buzsaki 2016) 
%                   -Required string-value pair inputs:
%                   'eventFiringRateMatrix', 'rateByPositionMatrix',
%                   'eventBinDuration'
%
% Inputs:
%       Positional: 
%           String indicating the shuffle type (see above)
%       String-value pair arguments: (different shuffle types require
%               a different subset of inputs, or all inputs can be
%               provided)
%           'bayesPostProb' - [nTimeBin X nSpatialBin] matrix of Posterior
%                   probabilities (see 'placeBayesLogBuffered')
%           'binIDs' - [nTimeBin X 1] column vector of the unique numeric
%                   ids indicating which event each of the time bins
%                   correspond to.
%           'eventFiringRateMatrix' - [nTimeBin X nCell] matrix of
%                   within-event by-time-bin firing rates
%           'rateByPositionMatrix' - [nCell X nSpatialBin] matrix of firing
%                   rates by position (place fields)
%           'eventBinDuration' - scalar duration (in seconds) of the bins
%                   within each event.
%
% Outputs:
%       'shuffledEvents' -  [nTimeBin X nSpatialBin] matrix of shuffled
%               Posterior probabilities
%       'newRateMatrix' - for shuffle types that create new or shuffled
%               event firing rates ('withinEventActiveCellIDPermutation &
%               'crossEventPoissonSurrogate') the shuffled rate matrix is
%               also output, in other cases this output is left empty.
% 
% Written by Andres Grosmark in 2021

bayesPostProb = [];
binIDs = [];
eventFiringRateMatrix = [];
rateByPositionMatrix = [];
eventBinDuration = [];
for i = 1:(length(varargin) - 1)
    if ischar(varargin{i})
        switch lower(varargin{i})
            case lower('bayesPostProb')
                bayesPostProb = varargin{i + 1};
            case lower('binIDs')
                binIDs = varargin{i + 1};
            case lower('eventFiringRateMatrix')
                eventFiringRateMatrix = varargin{i + 1};
            case lower('rateByPositionMatrix')
                rateByPositionMatrix = varargin{i + 1};
            case lower('eventBinDuration')
                eventBinDuration = varargin{i + 1};
        end
    end
end

newRateMatrix = [];
switch lower(shuffleType)
    case lower('timeBinPermutation')
        %check that we have the neccesary inputs:
        if isempty(bayesPostProb)
            error('In order to perform the ''timeBinPermutation'' a valid ''bayesPostProb'' matrix must be inputed as a value-pair argument');
        end
        if isempty(binIDs)
            error('In order to perform the ''timeBinPermutation'' a valid ''binIDs'' matrix must be inputed as a value-pair argument');
        end
        
        shuffledEvents = NaN(size(bayesPostProb));
        uEvents = unique(binIDs);
        for E = 1:length(uEvents)
            fE = find(binIDs == uEvents(E));
            shuffledEvents(fE, :) = bayesPostProb(fE(randperm(length(fE))), :);      % randomly permute bins within each event
        end
    case lower('withinEventActiveCellIDPermutation')
        %check that we have the neccesary inputs:
        if isempty(eventFiringRateMatrix)
            error('In order to perform the ''withinEventActiveCellIDPermutation'' a valid ''eventFiringRateMatrix'' matrix must be inputed as a value-pair argument');
        end
        if isempty(binIDs)
            error('In order to perform the ''withinEventActiveCellIDPermutation'' a valid ''binIDs'' matrix must be inputed as a value-pair argument');
        end
        if isempty(rateByPositionMatrix)
            error('In order to perform the ''withinEventActiveCellIDPermutation'' a valid ''rateByPositionMatrix'' matrix must be inputed as a value-pair argument');
        end
        if isempty(eventBinDuration)
            error('In order to perform the ''withinEventActiveCellIDPermutation'' a valid ''eventBinDuration'' matrix must be inputed as a value-pair argument');
        end
        
        uEvents = unique(binIDs);
        newRateMatrix = zeros(size(eventFiringRateMatrix));
        for i = 1:size(uEvents, 1)
            activeCells = find(sum(eventFiringRateMatrix(binIDs == uEvents(i), :), 1) > 0);
            newRateMatrix(binIDs == uEvents(i), activeCells) = eventFiringRateMatrix(binIDs == uEvents(i), activeCells(randperm(length(activeCells))));
        end
        shuffledEvents = placeBayesLogBuffered(newRateMatrix, rateByPositionMatrix, eventBinDuration);
    case lower('crossEventPoissonSurrogate')
        %check that we have the neccesary inputs:
        if isempty(eventFiringRateMatrix)
            error('In order to perform the ''crossEventPoissonSurrogate'' a valid ''eventFiringRateMatrix'' matrix must be inputed as a value-pair argument');
        end
        if isempty(rateByPositionMatrix)
            error('In order to perform the ''crossEventPoissonSurrogate'' a valid ''rateByPositionMatrix'' matrix must be inputed as a value-pair argument');
        end
        if isempty(eventBinDuration)
            error('In order to perform the ''crossEventPoissonSurrogate'' a valid ''eventBinDuration'' matrix must be inputed as a value-pair argument');
        end
        spikeCounts = eventFiringRateMatrix*eventBinDuration;
        newRateMatrix = NaN(size(eventFiringRateMatrix));
        for i = 1:size(newRateMatrix, 2)
            newRateMatrix(:, i) = poissrnd(mean(spikeCounts(:, i)), [size(eventFiringRateMatrix, 1), 1])/eventBinDuration;
        end
        shuffledEvents = placeBayesLogBuffered(newRateMatrix, rateByPositionMatrix, eventBinDuration);
    case lower('circularPlaceFieldShuffle')
        %check that we have the neccesary inputs:
        if isempty(eventFiringRateMatrix)
            error('In order to perform the ''circularPlaceFieldShuffle'' a valid ''eventFiringRateMatrix'' matrix must be inputed as a value-pair argument');
        end
        if isempty(rateByPositionMatrix)
            error('In order to perform the ''circularPlaceFieldShuffle'' a valid ''rateByPositionMatrix'' matrix must be inputed as a value-pair argument');
        end
        if isempty(eventBinDuration)
            error('In order to perform the ''circularPlaceFieldShuffle'' a valid ''eventBinDuration'' matrix must be inputed as a value-pair argument');
        end
       rateByPositionMatrixShuffle = NaN(size(rateByPositionMatrix));
       for i = 1:size(rateByPositionMatrixShuffle, 1)
           rateByPositionMatrixShuffle(i, :) = circshift(rateByPositionMatrix(i, :), randi(size(rateByPositionMatrix, 2)), 2);
       end
       shuffledEvents = placeBayesLogBuffered(eventFiringRateMatrix, rateByPositionMatrixShuffle, eventBinDuration);
end
