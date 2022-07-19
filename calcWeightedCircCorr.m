function [circCorrCoeff] = calcWeightedCircCorr(xy, w)
%function out = calcWeightedCircCorr(xy, w)
% Calculate the correlation between one linear variable (x) and one
% circular variable (y) as weighted by the weight matrix (w)
% Input: 
%   'xy': a two column vector where each row corresponds to the particular
%       [x, y] coordinate of one element in the linearized weight vector w
%       Note: y must be circular variable in radians (between 0 and 2*pi) 
%   'w': a column vector of linearized weights (typically the previously
%       computed Bayesian posterior probability estimates)
%
% Output:
%   'circCorrCoeff' = non-negative scalar, weighted circo-linear coefficent 
%
% See: Grosmark, A, Sparks, F.T., Davis, M.J. & Losonczy, A. Reactivation
%       predicts the consolidation of unbiased long-term cognitive maps.
%       Nature Neuroscience (2021).
% See also: 'circ_corrcl' from the Circular Statistics Toolbox for Matlab
%
% Note: the current method for estimating circo-linear correlations does
% not resolve the sign of this relationship. For this reason it is
% currently used in conjuction with the Line-Casting (Radon transform)
% approach. In the future perhaps the weighted correlation approach can be
% combined with the circo-linear method developed in the following paper
% which does resolve this sign. This would likely be more computationally
% intensive though:
%   Kempter R, Leibold C, Buzsáki G, Diba K, Schmidt R. Quantifying
%   circular-linear associations: hippocampal phase precession. J Neurosci
%   Methods. 2012 May 30;207(1):113-24. Epub 2012 Apr 7.
%
% Written by Andres Grosmark in 2021


rxs = calcWeightedLinearCorr([xy(:, 1), sin(xy(:, 2))], w);
rxc = calcWeightedLinearCorr([xy(:, 1), cos(xy(:, 2))], w);
rcs = calcWeightedLinearCorr([sin(xy(:, 2)),cos(xy(:, 2))], w);


circCorrCoeff = sqrt((rxc^2 + rxs^2 - 2*rxc*rxs*rcs)/(1-rcs^2)); 


