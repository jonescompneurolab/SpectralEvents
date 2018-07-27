% This test script implements the main pipeline of the spectralevents 
% toolbox. All data is taken from sessions used in Shin et al. eLife 2017
% (http://dx.doi.org/10.5061/dryad.pn931).
% Ryan Thorpe, Human Electrophysiology and Computational Neuroscience Lab, 
% Brown University

clear all

% Load data sessions from the same experimental setup so that spectral 
% event features are characterized for a homogenous set of experimental
% conditions
numSubj = 10;
x = cell(1,numSubj);
classLabels = cell(1,numSubj);
for subj_i=1:numSubj
    load(['test_data/prestim_humandetection_600hzMEG_subject',num2str(subj_i),'.mat'])
    x{subj_i} = prestim_raw_yes_no'; %Raw timeseries for detection/non-detection prestimulus LFP
    classLabels{subj_i} = YorN'; %Vector of trial classification labels
    clear prestim_raw_yes_no YorN
end

% Set dataset and analysis parameters
eventBand = [15,29]; %Frequency range of spectral events
vis = true; %Generate standard visualization plots for event features across all subjects/sessions
fVec = (1:60); %Vector of fequency values over which to calculate TFR
Fs = 600; %Sampling rate of time-series
%tVec = (0:1/Fs:1);
[specEvents,TFRs,X] = spectralevents(eventBand,vis,fVec,Fs,x,classLabels); %Run spectral event analysis

