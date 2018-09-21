% This script tests and demonstrates implementation of the main pipeline of
% the spectralevents toolbox. The datasets referenced below are from 
% sessions used in Shin et al. eLife 2017 
% (http://dx.doi.org/10.5061/dryad.pn931).

%   -----------------------------------------------------------------------
%   SpectralEvents::test
%   Copyright (C) 2018  Ryan Thorpe
%
%   This file is part of the SpectralEvents toolbox.
% 
%   SpectralEvents is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   SpectralEvents is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.
%   -----------------------------------------------------------------------

clear all, close all

% Load data sessions/subjects from the same experimental setup so that 
% spectral event features are differentially characterized only between the 
% desired trial classification labels: in this case, detection vs. 
% non-detection
numSubj = 10;
x = cell(1,numSubj);
classLabels = cell(1,numSubj);
for subj_i=1:numSubj
    load(['test_data/prestim_humandetection_600hzMEG_subject',num2str(subj_i),'.mat'])
    x{subj_i} = prestim_raw_yes_no'; %Time-by-trial matrix of timeseries trials for detection/non-detection prestimulus MEG
    classLabels{subj_i} = YorN'; %Column vector of trial classification labels
    clear prestim_raw_yes_no YorN
end

% Set dataset and analysis parameters
eventBand = [15,29]; %Frequency range of spectral events
fVec = (1:60); %Vector of fequency values over which to calculate TFR
Fs = 600; %Sampling rate of time-series
findMethod = 1; %Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
vis = true; %Generate standard visualization plots for event features across all subjects/sessions
%tVec = (1/Fs:1/Fs:1);
[specEvents,TFRs,timeseries] = spectralevents(eventBand,fVec,Fs,findMethod,vis,x,classLabels); %Run spectral event analysis
