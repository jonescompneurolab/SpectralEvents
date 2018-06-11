clear all
load('/Users/rv_01/Documents/research/data/data_shin_et_al_2018/prestim_mousedetection_mouse1_session2.mat')

thrFOM=6; %Threshold Factors of Median
rhythm_band_inds = 15:29; % indices of fVec that corresponds to 15 to 29 Hz band.
%rhythmid='beta';

s = SpectralEventAnalysis(rhythm_band_inds, thrFOM, tVec, fVec, prestim_TFR_yes_no, YorN);