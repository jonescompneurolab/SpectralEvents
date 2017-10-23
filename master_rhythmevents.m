resultPath=uigetdir('', 'folder to load data from, and save data to'); % user can input where 
thrFOM=6;
rhythm_band_inds = 15:29; % indices of fVec that corresponds to 15 to 29 Hz band.
rhythmid='beta';

% dataset type : 1 for mouse detection, 2 for human detection, 3 for human attention
%% mouse preprocess TFR
datatypeid='mousedetection';
subject_list={'mouse1_session1', 'mouse1_session2', 'mouse1_session3', 'mouse1_session4', 'mouse1_session5', ....
  'mouse2_session1', 'mouse2_session2', 'mouse2_session3', 'mouse2_session4', 'mouse2_session5'};

for s=1:10
  clearvars -except subject_list subjind thrFOM rhythm_band_inds rhythmid  
  
  load(strcat(resultPath, '\prestim_', datatypeid, '_', subject_list{s}, '.mat'))
    
  indsoi = 1:length(tVec); % Indices correspond to t = -1000:0 ms.
  rhythm_localmax_analysis(resultPath, datatypeid, subject_list, s, rhythmid, rhythm_band_inds, tVec, fVec, indsoi, prestim_TFR_yes_no, YorN)
  rhythm_event_analysis(rhythmid, datatypeid, resultPath, subject_list, s, thrFOM, tVec, fVec, prestim_TFR_yes_no, YorN)
end

%% human detection data
datatypeid='humandetection';
subject_list = {'subject1', 'subject2', 'subject3', 'subject4', 'subject5', 'subject6', 'subject7', 'subject8', 'subject9', 'subject10'}; % List of Subject ID Numbers

for s=1:10
  clearvars -except subject_list subjind thrFOM rhythm_band_inds rhythmid  
  
  load(strcat(resultPath, '\prestim_', datatypeid, '_', subject_list{s}, '.mat'))  
    
  indsoi = 1:length(tVec); % Indices correspond to t = -1000:0 ms.
  rhythm_localmax_analysis(resultPath, datatypeid, subject_list, s, rhythmid, rhythm_band_inds, tVec, fVec, indsoi, prestim_TFR_yes_no, YorN)  
  rhythm_event_analysis(rhythmid, datatypeid, resultPath, subject_list, s, thrFOM, tVec, fVec, prestim_TFR_yes_no, YorN)
end

%% human attention data: detect trials only
datatypeid='humanattention';
subject_list = {'111' '112' '113' '115' '117' '118' '119' '120' '121' '122'}; % List of Subject ID Numbers

for s=1:10
  clearvars -except subject_list subjind thrFOM rhythm_band_inds rhythmid  
  
  load(strcat(resultPath, '\prestim_', datatypeid, '_', subject_list{s}, '.mat'))
     
  indsoi = 1:length(tVec); % Indices correspond to t = -1000:0 ms.
  rhythm_localmax_analysis(resultPath, datatypeid, subject_list, s, rhythmid, rhythm_band_inds, tVec, fVec, indsoi, prestim_TFR_yes_no, YorN)
  rhythm_event_analysis(rhythmid, datatypeid, resultPath, subject_list, s, thrFOM, tVec, fVec, prestim_TFR_yes_no, YorN)
end
