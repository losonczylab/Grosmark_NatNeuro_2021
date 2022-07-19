# Linear and Circular Trajectory Replay Analysis

Author: Andres Grosmark 2021

These functions were used to assess sequence replay in the study:

*Grosmark AD, Sparks FT, Davis MJ, Losonczy A. Reactivation predicts the consolidation of unbiased long-term cognitive maps. Nat Neurosci. 2021 Nov;24(11):1574-1585.*

These Matlab functions perform:
1) Bayesian decoding analysis (Zhang et al. 1998) 
2) Linear Weighted correlation replay analysis (Wu & Foster, 2014), 
3) Circo-Linear weighted correlation replay analysis (Grosmark et al. 2021)
4) Line casting (or Radon) replay analysis (see Davidson et al. 2009, here also adapted to circular replay).

In Grosmark 2021 these methods were applied to Ca2+ data after spike-deconvolution and thresholding (see Methods) and are designed to be equally applicable spike-deconvolved Ca2+ analysis as well as to electrophysiologically recorded spike analysis - they are not directly suited for the analysis of continous signals (for instance, Df over F).

## Requirements:

Matlab (tested on versions 2018a and 2022a), and the Matlab ‘Image Processing’ and ‘Statistics and Machine Learning’ toolboxes.

## How to run this code:

Clone or download this repository then add the folder to your path (or run files from within folder). The functionality of these analysis functions is demonstrated by two ‘demo’ scripts:

1.  **‘demo_LinearReplayAnalysis.m’** – performs linear (weighted and line-casting) replay analysis on 5 synthetic replay event created by ‘createSyntheticEvents.m’ and plots the results.
2.  *‘***demo_CircularReplayAnalysis.m***’* – performs linear (weighted and line-casting) replay analysis on 5 synthetic replay event created by ‘createSyntheticEvents.m’ and plots the results.

## Key functions:

-**'placeBayesLogBuffered’** – performs Bayesian decoding analysis

-**'calcWeightedLinearCorr’** – performs linear weighted correlation analysis

-**'calcWeightedCircCorr’** - performs circo-linear weighted correlation analysis

-**'calcRadonReplay’** – performs ‘line-casting’ (Radon transformation) replay analysis

Note that unlike the linear correlation analysis, the results of the circo-linear correlation analysis are not signed (i.e. are always positive) and therefore the Radon transform analysis is used to disambiguate between forward and reverse replay events.

Input-output details, notes and citations are provided in the help for each of the functions.
