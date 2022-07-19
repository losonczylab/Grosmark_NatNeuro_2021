function out = calcWeightedLinearCorr(xy, w);
% function out = calcWeightedLinearCorr(xy, w);
% Calculate linear weighted correlations
% Inputs:
%   'xy': a two column vector where each row corresponds to the particular
%       [x, y] coordinate of one element in the linearized weight vector w
%   'w': a column vector of linearized weights (typically the previously
%       computed Bayesian posterior probability estimates)
% Output: 
%   'weighted correlation coefficent' (scalar).
%
%  See: Zhang K, Ginzburg I, McNaughton BL, Sejnowski TJ. Interpreting
%       neuronal population activity by reconstruction: unified framework
%       with application to hippocampal place cells. J Neurophysiol. 1998
%       Feb;79(2):1017-44. For original implementation.
%
% Written by Andres Grosmark in 2021

mxy = sum(xy.*repmat(w, 1, 2)/sum(w));
covxy = sum(w.*(xy(:, 1) - mxy(1)).*(xy(:, 2) - mxy(2)))/sum(w);
covxx = sum(w.*(xy(:, 1) - mxy(1)).*(xy(:, 1) - mxy(1)))/sum(w);
covyy = sum(w.*(xy(:, 2) - mxy(2)).*(xy(:, 2) - mxy(2)))/sum(w);
out = covxy/(sqrt(covyy*covxx));