% This script tests and demonstrates implementation of the main pipeline of
% the spectralevents toolbox. The datasets referenced below are from 
% sessions used in Shin et al. eLife 2017 
% (http://dx.doi.org/10.5061/dryad.pn931).

clear all, close all

% Load data sessions/subjects from the same experimental setup so that 
% spectral event features are differentially characterized only between the 
% desired trial classification labels: in this case, detection vs. 
% non-detection
numSubj = 10;
x = cell(1,numSubj);
classLabels = cell(1,numSubj);
for subj_i=1:numSubj
    load(['data/prestim_humandetection_600hzMEG_subject',num2str(subj_i),'.mat'])
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

% Save figures
classes = [0,1];
for subj_i=1:numSubj
    for class_i=1:2
        figName = strcat('./example_output/matlab/prestim_humandetection_600hzMEG_subject', num2str(subj_i), '_class_', num2str(classes(class_i)), '.png');
        saveas(figure((subj_i-1)*2+class_i),figName);
    end
end