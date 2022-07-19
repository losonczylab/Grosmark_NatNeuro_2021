function [shuffledEvents, newRateMatrix] = shufflePopulationEvents(shuffleType, varargin)
%% shuffleTypes = 'timeBinPermutation', 'withinEventActiveCellIDPermutation', 'crossEventPoissonSurrogate', 'circularPlaceFieldShuffle'

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






