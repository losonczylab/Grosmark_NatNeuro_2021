%% Demo script for linear replay analysis
% Written by Andres Grosmark in 2021

rng(1); % for consistancy.

totalMazeLength = 2; %total circumferance of ciruclar maze in meters

nShuffles = 500; % number of shuffles to perform

radonBoxCarBinWidth = 5; % size of box-car smoothing window applied to 
% data before radon analysis (see below)


%% Note: comment/uncomment to select shuffle type (out of the 4
% implemented options, see the 'shufflePopulationEvents' function for more
% details).

shuffleType = 'timeBinPermutation';
% shuffleType = 'withinEventActiveCellIDPermutation';
% shuffleType = 'crossEventPoissonSurrogate';
% shuffleType = 'circularPlaceFieldShuffle';

%% Create the synthetic events to test out the replay analysis

synthEvents = createSyntheticEvents('nTimeBinsPerEvent', 15);
LinReplayOutput = [];

%% Plot synthetic place fields

figure; 
cellID = 1:size(synthEvents.rateByPositionMatrix, 1);
spatialBinCenters = linspace(0, totalMazeLength, size(synthEvents.rateByPositionMatrix, 2) + 1);
spatialBinCenters = spatialBinCenters(1:(end - 1)) + (totalMazeLength/(2*size(synthEvents.rateByPositionMatrix, 2)));
imagesc(spatialBinCenters, cellID, synthEvents.rateByPositionMatrix)
cBar = colorbar;
cBar.Label.String = 'Firing Rate (Hz)';
cBar.Label.Rotation = -90;
cBar.Label.VerticalAlignment = 'bottom';
xlabel('Position (m)');
ylabel('Cell #');
title('Synthetic Place Fields Map')

%% Calculate the Bayesian posteriors for each bin in the synthetic events using the the synthetic  
% place fields as a template
LinReplayOutput.bayesPostProb = placeBayesLogBuffered(synthEvents.eventFiringRateMatrix, ...
    synthEvents.rateByPositionMatrix, synthEvents.params.eventBinDuration);

LinReplayOutput.binIDs = synthEvents.eventFiringRateMatrixID;

%% Create helper matrix to save time with weighted correlation step below
% xyByBin is an [nx2] matrix of the linearized x (space) and y (time) positions
% for each element of an event in 'LinReplayOutput.bayesPostProb'.
% Because (assuming the all the events have the same number of total
% position bins) this x-y matrix is shared across events with the same
% number of time bins, we pre-calculate it for each unique number of
% time-bins. This saves a lot of unnecessary repetition of this calculation
% down the road, and thus saves computational time. xyByBin is neccesary
% for running the weighted correlation functions which calculate the actual
% replay scores.

h = histc(synthEvents.eventFiringRateMatrixID, unique(synthEvents.eventFiringRateMatrixID));
uBins = unique(h);

ySpatial = linspace(0, totalMazeLength,  synthEvents.params.nSpatialBins + 1);
ySpatial = ySpatial(1:synthEvents.params.nSpatialBins);
xyByBin = {};
for i = 1:length(uBins)
    x = repmat((1:uBins(i))', [1, 100]);
    y = repmat(ySpatial, [uBins(i), 1]);
    xyByBin{uBins(i)} = [reshape(x, [], 1), reshape(y, [], 1)];
end

%% Calculate linear weighted correlations
uEvents = unique(synthEvents.eventFiringRateMatrixID);
LinReplayOutput.eventIDs = uEvents;
LinReplayOutput.linearWeightedCorrCoeff = [];
for E = 1:length(uEvents)
    bayesPost = LinReplayOutput.bayesPostProb(synthEvents.eventFiringRateMatrixID == E, :);
    nBins = sum(synthEvents.eventFiringRateMatrixID == uEvents(E));
    linearCorr = calcWeightedLinearCorr(xyByBin{nBins}, bayesPost(:));
    LinReplayOutput.linearWeightedCorrCoeff = [LinReplayOutput.linearWeightedCorrCoeff; linearCorr];
end

%% Calculate maximum jump distance in each event
uEvents = unique(synthEvents.eventFiringRateMatrixID);
spatialBinLength = totalMazeLength./synthEvents.params.nSpatialBins;
LinReplayOutput.MaxJumpDist = [];
for E = 1:length(uEvents)
    bayesPost = LinReplayOutput.bayesPostProb(synthEvents.eventFiringRateMatrixID == E, :);
    [~, maxPosBin] = max(bayesPost, [], 2);
    maxJump = max(abs(diff(maxPosBin)))*spatialBinLength;
    LinReplayOutput.MaxJumpDist = [LinReplayOutput.MaxJumpDist; maxJump];
end

%% Pre-compute a look-up table containing the information of all the line 'cast' by the
% Radon transform (used in Radon-replay calculation) to save time later on
% (this makes a bigger difference when lots of events are being analyzed
% and less for these synthetic examples). 

radonLookupTable = makeRadonLookupTable(size(LinReplayOutput.bayesPostProb, 2), uBins);

%% Calculate radon ('line-casting') replay - this is useful for estimating the exact replay trajectory
% Note: this method is applied to a (within-temporal bin) smoothed version
% of the posterior probability matrix, allowing for slight non-linearities
% in the replay signal recoverable with this method. Note that this isn't
% neccesary for the current, synthetic, set of examples but can be useful
% for real data.
if isempty(radonBoxCarBinWidth)    
    [LinReplayOutput.Radon.posMean, LinReplayOutput.Radon.id, LinReplayOutput.Radon.slope, ...
        LinReplayOutput.Radon.pathLength, LinReplayOutput.Radon.pointsXYXY] = ...
        calcRadonReplay( LinReplayOutput.bayesPostProb', synthEvents.eventFiringRateMatrixID, radonLookupTable);
else
    LinReplayOutput.Radon = [];
    boxCar1 = ones(radonBoxCarBinWidth, 1);
    bayesReconSmoothed = convWith(LinReplayOutput.bayesPostProb', boxCar1);
    bayesReconSmoothed = bayesReconSmoothed./repmat(sum(bayesReconSmoothed, 1), [size(bayesReconSmoothed, 1), 1]);
    
    [LinReplayOutput.Radon.posMean, LinReplayOutput.Radon.id, LinReplayOutput.Radon.slope,...
        LinReplayOutput.Radon.pathLength, LinReplayOutput.Radon.pointsXYXY] = ...
        calcRadonReplay(bayesReconSmoothed, synthEvents.eventFiringRateMatrixID, radonLookupTable);
end

%% convert the radon best fit line and slope output units from bins to meters and seconds:

LinReplayOutput.Radon.pointsXYXY(:, [2, 4]) = (LinReplayOutput.Radon.pointsXYXY(:, [2, 4]) - 0.5)...
    *(totalMazeLength/synthEvents.params.nSpatialBins); % convert from spatial bins to distance

LinReplayOutput.Radon.pointsXYXY(:, [1, 3]) = (LinReplayOutput.Radon.pointsXYXY(:, [1, 3]) - 0.5)...
    *synthEvents.params.eventBinDuration; % convert from temporal bins to time

LinReplayOutput.Radon.slope = LinReplayOutput.Radon.slope...
    *(totalMazeLength/synthEvents.params.nSpatialBins)/synthEvents.params.eventBinDuration;


%% Calculate shuffled weighted correlations and maximum jump distances for each event

uEvents = unique(synthEvents.eventFiringRateMatrixID);

LinReplayOutput.Shuffle = [];
LinReplayOutput.Shuffle.linearWeightedCorrCoeff = NaN(length(uEvents), nShuffles);
LinReplayOutput.Shuffle.maxJump = NaN(length(uEvents), nShuffles);
for sh = 1:nShuffles
    switch lower(shuffleType) %different types of shuffle require different inputs (one can also just pass alll the inputs instead)
        case lower('timeBinPermutation')
            shuffledBayesianPosteriors = shufflePopulationEvents(shuffleType, 'bayesPostProb', ...
                LinReplayOutput.bayesPostProb, 'binIDs', LinReplayOutput.binIDs);
        case lower('withinEventActiveCellIDPermutation')
            shuffledBayesianPosteriors = shufflePopulationEvents(shuffleType, 'eventFiringRateMatrix', ...
                synthEvents.eventFiringRateMatrix, 'binIDs', synthEvents.eventFiringRateMatrixID, ...
                'rateByPositionMatrix', synthEvents.rateByPositionMatrix, 'eventBinDuration', ...
                synthEvents.params.eventBinDuration);
        case lower('crossEventPoissonSurrogate')
            shuffledBayesianPosteriors = shufflePopulationEvents(shuffleType, 'eventFiringRateMatrix', ...
                synthEvents.eventFiringRateMatrix, 'rateByPositionMatrix', synthEvents.rateByPositionMatrix, ...linearcalcWeightedLinearCorr
                'eventBinDuration', synthEvents.params.eventBinDuration);
        case lower('circularPlaceFieldShuffle')
            shuffledBayesianPosteriors = shufflePopulationEvents(shuffleType, 'eventFiringRateMatrix', ...
                synthEvents.eventFiringRateMatrix, 'rateByPositionMatrix', synthEvents.rateByPositionMatrix, ...
                'eventBinDuration', synthEvents.params.eventBinDuration);
    end
    
    for E = 1:length(uEvents)
        nBins = sum(synthEvents.eventFiringRateMatrixID == uEvents(E));
        bayesPost = shuffledBayesianPosteriors(synthEvents.eventFiringRateMatrixID == uEvents(E), :);
        linearCorr = calcWeightedLinearCorr(xyByBin{nBins}, bayesPost(:));
        LinReplayOutput.Shuffle.linearWeightedCorrCoeff(E, sh) = linearCorr;
        
        [~, maxPosBin] = max(bayesPost, [], 2);
        maxJump = max(abs(diff(maxPosBin)))*spatialBinLength;
        LinReplayOutput.Shuffle.maxJump(E, sh) = maxJump;
    end
end

%% calculate weighted correltation p-values and 'rZ' scores
LinReplayOutput.linearWeightedCorr_PVal = mean(abs(LinReplayOutput.Shuffle.linearWeightedCorrCoeff) ...
    >= abs(repmat(LinReplayOutput.linearWeightedCorrCoeff, [1, nShuffles])), 2);
LinReplayOutput.linearWeightedCorr_rZ = (abs(LinReplayOutput.linearWeightedCorrCoeff) ...
    - mean(abs(LinReplayOutput.Shuffle.linearWeightedCorrCoeff), 2))./...
    std(abs(LinReplayOutput.Shuffle.linearWeightedCorrCoeff), [], 2);

%% Print the linear replay analysis results
uEvents = unique(synthEvents.eventFiringRateMatrixID);
for E = 1:length(uEvents)
    switch sign(LinReplayOutput.linearWeightedCorrCoeff(E))
        case 1
            repType = 'Forward';
        case -1
            repType = 'Reverse';
    end
    display(['Event #', int2str(uEvents(E)), ', Synthetic Event Type: ', synthEvents.eventType{uEvents(E)}]);
    display(['Linear weighted correlation coefficent: ', num2str(LinReplayOutput.linearWeightedCorrCoeff(E), 5)]);
    display(['Event p = ', num2str(LinReplayOutput.linearWeightedCorr_PVal(E), 3), ', rZ = ', ...
        num2str(LinReplayOutput.linearWeightedCorr_rZ(E), 3)]);
    display(['Radon (line-casting) slope: ', num2str(LinReplayOutput.Radon.slope(E), 3), ...
        ' m/s, Putative ', repType, ' replay event.']);
    display(' ');
end

%% Plot the linear replay analysis results
spatialBinLength = totalMazeLength./synthEvents.params.nSpatialBins;
spatialBinCenters = linspace(0, totalMazeLength - spatialBinLength, synthEvents.params.nSpatialBins) ...
    + spatialBinLength/2;

eventDuration = synthEvents.params.nTimeBinsPerEvent*synthEvents.params.eventBinDuration;
timeBinCenter = linspace(0, eventDuration - synthEvents.params.eventBinDuration, ...
    synthEvents.params.nTimeBinsPerEvent) + synthEvents.params.eventBinDuration/2;
cellId = 1:synthEvents.params.nPlaceCells;

figure('Units', 'normalized', 'Position', [0.033854     0.062037      0.84583      0.81296], ...
    'Color', [0.94         0.94         0.94]);

xAxesPos = linspace(0.08, 0.82, 5);
allAxes = {};
lineOverlayColor = [0, 0.866, 0.0039]; % a nice darkish green.
for E = 1:5
    % first row of axes: per event firing rates
    allAxes{1, E} = axes(); 
    imagesc(timeBinCenter, cellId,  ...
        synthEvents.eventFiringRateMatrix(synthEvents.eventFiringRateMatrixID == uEvents(E), :));
    cBar = colorbar;
    cBar.Label.String = 'Within-Bin Firing Rate (Hz)';
    cBar.Label.Rotation = -90;
    cBar.Label.VerticalAlignment = 'bottom';
    formattedTitle = {};
    for i = 1:length(synthEvents.eventTypeMultiLine{E})
        formattedTitle{i} = ['\fontsize{12}\bf', synthEvents.eventTypeMultiLine{E}{i}];
    end
    title(formattedTitle);
    if E == 1
        ylabel('Cell #');
    end
    xlabel('Time (sec)');
    allAxes{1, E}.Position(1) = xAxesPos(E);
    allAxes{1, E}.Position(2) = 0.66;
    allAxes{1, E}.Position(3) = 0.11;
    allAxes{1, E}.Position(4) = 0.216;
    
    
    % second row of axes: per event Bayesian reconstructions
    allAxes{2, E} = axes(); 
    bayesPost = LinReplayOutput.bayesPostProb(synthEvents.eventFiringRateMatrixID == uEvents(E), :)';
    imagesc(timeBinCenter, spatialBinCenters, bayesPost);
    
    colormap(allAxes{2, E}, 'hot');
    cBar = colorbar;
    cBar.Label.String = 'Bayesian Post. Prob.';
    cBar.Label.Rotation = -90;
    cBar.Label.VerticalAlignment = 'bottom';
    hold on;
    radonPointsXYXY = LinReplayOutput.Radon.pointsXYXY(E, :);

    plot(radonPointsXYXY([1, 3]), radonPointsXYXY([2, 4]), ...
        'Color', lineOverlayColor, 'LineWidth', 1);
        
    if E == 1
        ylabel('Spatial Location (m)');
    end
    xlabel('Time (sec)');
    switch sign(LinReplayOutput.linearWeightedCorrCoeff(E))
        case 1
            repType = 'Forward';
        case -1
            repType = 'Reverse';
    end
        
    
    title({['\rm\fontsize{9}Linear Weighted Corr = ', num2str(LinReplayOutput.linearWeightedCorrCoeff(E), 3)], ...
        ['Max Jump. Dist = ', num2str(LinReplayOutput.MaxJumpDist(E), 3), ' m'], ...
        ['Radon Slope: ', num2str(LinReplayOutput.Radon.slope(E), 3), ' m/s, (', repType, ' rep.)']});
    
    allAxes{2, E}.Position(1) = xAxesPos(E);
    allAxes{2, E}.Position(2) = 0.3125;
    allAxes{2, E}.Position(3) = 0.11;
    allAxes{2, E}.Position(4) = 0.216;
    
    % third row second row of axes: per event observed vs. shuffled weighted correlation coefficents
    allAxes{3, E} = axes(); 
    allAxes{3, E}.Position(1) = xAxesPos(E);
    allAxes{3, E}.Position(2) = 0.06;
    allAxes{3, E}.Position(3) = 0.11;
    allAxes{3, E}.Position(4) = 0.15;
    hist(abs(LinReplayOutput.Shuffle.linearWeightedCorrCoeff(E, :)), 20);
    yL = get(gca, 'YLim');
    hold on; plot(repmat(abs(LinReplayOutput.linearWeightedCorrCoeff(E)), [1, 2]), yL, 'm', 'LineWidth', 1.5);
    ylim(yL);
    ylabel('# of shuffled events');
    xlabel('Weighted Corr. Coeff.');
    if LinReplayOutput.linearWeightedCorr_PVal(E) < 0.05
        title({['\bfEvent p-value = ', num2str(LinReplayOutput.linearWeightedCorr_PVal(E), 3)], ...
            ['\rmrZ = ', num2str(LinReplayOutput.linearWeightedCorr_rZ(E), 3)]});
    else
        title({['\rmEvent p-value = ', num2str(LinReplayOutput.linearWeightedCorr_PVal(E), 3)], ...
            ['\rmrZ = ', num2str(LinReplayOutput.linearWeightedCorr_rZ(E), 3)]});
    end    
end

% Title the figure:
TitleBox = uicontrol('style','text');
TitleBox.String = 'Linear Replay Analysis';
TitleBox.Units = 'normalized';
% TitleBox.BackgroundColor = 'None';
TitleBox.HorizontalAlignment = 'left';
TitleBox.FontAngle = 'italic';
TitleBox.FontWeight = 'bold';
TitleBox.FontSize = 16;
TitleBox.Position(1) = 0.07;
TitleBox.Position(2) = 0.935;
TitleBox.Position(3) = 0.2;
TitleBox.Position(4) = 0.05;
