function [Pr, prMax] = placeBayesLogBuffered1(Cr, rateMap, binLength)
%function [Pr, prMax] = placeBayesLogBuffered1(Cr, rateMap, binLength)
%Inputs: Cr - an [nTime-bin X nCell] rate matrix, Bayesian decoding is 
%             performed on each time bin
%        rateMap - [nCell x nSpatial-bin] rate matrix - this is the 
%             'template' used for the decoding
%        binLength - the length (in seconds) of each Time-bin in Cr
%Outputs: Pr - [nTime-bin X nSpatial-bin] Posterior probability matrix
%         Pr-Max = max(Pr, [], 2);
%
% see:
%   T. J. Davidson, F. Kloosterman, M. A. Wilson, Hippocampal replay of extended experience.
%       Neuron 63, 497–507 (2009). for original implementation. 
% 
% Some of the calculations are performed in log-space (and then
% exponantiated) to prevent numerical overflow issues
% The analysis is performed in a 'chunk-wise' manner which I have found to
% be slightly faster (the size of the chunks are set by the 'buffer'
% variable).
%
% Written by Andres Grosmark in 2021
buffer = 12;

if size(Cr, 2) ~= size(rateMap, 1)
    error('Cr must be a [nTime-bin x nCell] matrix and rateMap must be [nCell x nSpatial-bins]');
end

Cr = Cr*binLength;
rateMap = rateMap' + (10^-10); % this small offset is neccesary to prevent -Inf results 
% when we take the log of this matrix below

term2 = ((-1)*binLength*sum(rateMap'));
Pr = [];
Pr(size(Cr, 1), size(rateMap, 1)) = NaN;
prMax = [];
prMax(size(Cr, 1), 1) = 0;

% These permutations are some tricks to make the calculations below a bit
% faster.
b = permute(log(rateMap'), [3, 2, 1]);
Cr = permute(Cr, [1, 3, 2]);

for S = 1:buffer:size(Cr, 1)
    T = S + buffer - 1;
    if T > size(Cr, 1)
        T = size(Cr, 1);
    end
    
    c =  repmat(Cr(S:T, :, :), [1, size(rateMap, 1), 1]);  
    u = sum(c.*repmat(b, [size(c, 1), 1, 1]), 3) + repmat(term2, size(c, 1), 1);
    u = u - repmat(max(u, [], 2), 1, size(u, 2));
    
    Pr(S:T, :) = exp(u);  
    Pr(S:T, :) = Pr(S:T, :)./repmat(sum(Pr(S:T, :), 2), 1, size(Pr(S:T, :), 2));
    
end

[~, prMax] = nanmax(Pr, [], 2);
if sum(sum(isinf(Pr))) > 0
    error('Do Not Approach the Infitnite');
end

if sum(sum(isnan(Pr))) > 0
    error('What is ''not a nubmer''?');
end