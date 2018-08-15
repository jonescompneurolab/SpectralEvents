% This test script implements the main pipeline of the spectralevents 
% toolbox. The datasets referenced below are from sessions used in Shin et 
% al. eLife 2017 (http://dx.doi.org/10.5061/dryad.pn931).

clear all

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
vis = true; %Generate standard visualization plots for event features across all subjects/sessions
fVec = (1:60); %Vector of fequency values over which to calculate TFR
Fs = 600; %Sampling rate of time-series
%tVec = (1/Fs:1/Fs:1);
[specEvents,TFRs,timeseries] = spectralevents(eventBand,fVec,Fs,vis,x,classLabels); %Run spectral event analysis
