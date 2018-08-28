function specEv_struct = spectralevents_find_1(eventBand, thrFOM, tVec, fVec, TFR, classLabels)
% SPECTRALEVENTS_FIND_1 Algorithm for finding and calculating spectral 
%   events on a trial-by-trial basis of of a single subject/session. Uses 
%   the following method (used as the primary event detection method in 
%   Shin et al. eLife 2017):
%   1) Retrieve all local maxima in TFR using imregionalmax
%   2) Pick out maxima above threshold and within the frequency band of 
%      interest
%   3) Identify and organize event features
%
% specEv_struct = spectralevents_find_1(eventBand,thrFOM,tVec,fVec,TFR,classLabels)
% 
% Inputs:
%   eventBand - range of frequencies ([Fmin_event Fmax_event]; Hz) over 
%       which above-threshold spectral power events are classified.
%   thrFOM - factors of median threshold; positive real number used to
%       threshold local maxima and classify events (see Shin et al. eLife 
%       2017 for discussion concerning this value).
%   tVec - time vector (s) over which the time-frequency response (TFR) is 
%       calcuated.
%   fVec - frequency vector (Hz) over which the time-frequency response 
%       (TFR) is calcuated.
%   TFR - time-frequency response (TFR) (frequency-by-time-trial) for a
%       single subject/session.
%   classLabels - numeric or logical 1-row array of trial classification 
%       labels; associates each trial of the given subject/session to an 
%       experimental condition/outcome/state (e.g., hit or miss, detect or 
%       non-detect, attend-to or attend away).
%
% Outputs:
%   specEv_struct - event feature structure with three main sub-structures:
%   TrialSummary (trial-level features), Events (individual event 
%   characteristics), and IEI (inter-event intervals from all trials and 
%   those associated with only a given class label).
%
% See also SPECTRALEVENTS, SPECTRALEVENTS_FIND_2, SPECTRALEVENTS_VIS.

% Initialize general data parameters
eventBand_inds = fVec(fVec>=eventBand(1) & fVec<=eventBand(2)); %Indices of freq vector within eventBand
flength = size(TFR,1); %Number of elements in discrete frequency spectrum
tlength = size(TFR,2); %Number of points in time
numTrials = size(TFR,3); %Number of trials
classes = unique(classLabels);

% Validate consistency of parameter dimensions
if flength~=length(fVec) || tlength~=length(tVec) || numTrials~=length(classLabels)
  error('Mismatch in input parameter dimensions!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1) Retrieve all local maxima in TFR using imregionalmax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TFRlocalmax: 12 column matrix for storing local max event metrics: trial 
% index, hit/miss, maxima frequency, lowerbound frequency, upperbound 
% frequency, frequency span, maxima timing, event onset timing, event 
% offset timing, event duration, maxima power, maxima/median power
TFRlocalmax = [];

% Finds_localmax: stores peak frequecy at each local max (columns) for each
% trial (rows)
Finds_localmax = [];

% Retrieve all local maxima in TFR using imregionalmax
for ti=1:numTrials
  [peakF,peakT] = find(imregionalmax(squeeze(TFR(:,:,ti)))); %Indices of max local power
  peakpower = TFR(find(imregionalmax(squeeze(TFR(:,:,ti))))+(ti-1)*flength*tlength); %Power values at local maxima (vector; compiles across frequencies and time)
  
  % Find local maxima lowerbound, upperbound, and full width at half max
  % for both frequency and time
  Ffwhm = NaN(numel(peakpower),3); %2D matrix for freq-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
  Tfwhm = NaN(numel(peakpower),3); %2D matrix for time-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
  for lmi=1:numel(peakpower)
    lmF_underthr = find(squeeze(TFR(:,peakT(lmi),ti) < peakpower(lmi)/2)); %Indices of TFR frequencies of < half max power at the time of a given local peak
    if ~isempty(find(lmF_underthr < peakF(lmi), 1)) && ~isempty(find(lmF_underthr > peakF(lmi), 1))
      Ffwhm(lmi,1) = fVec(lmF_underthr(find(lmF_underthr < peakF(lmi),1,'last'))+1);
      Ffwhm(lmi,2) = fVec(lmF_underthr(find(lmF_underthr > peakF(lmi),1,'first'))-1);
      Ffwhm(lmi,3) = Ffwhm(lmi,2)-Ffwhm(lmi,1)+ min(diff(fVec));
    elseif isempty(find(lmF_underthr < peakF(lmi),1)) && ~isempty(find(lmF_underthr > peakF(lmi),1))
      Ffwhm(lmi,1) = fVec(1);
      Ffwhm(lmi,2) = fVec(lmF_underthr(find(lmF_underthr > peakF(lmi),1,'first'))-1);
      Ffwhm(lmi,3) = 2*(Ffwhm(lmi,2)-fVec(peakF(lmi)))+ min(diff(fVec));
    elseif ~isempty(find(lmF_underthr < peakF(lmi),1)) && isempty(find(lmF_underthr > peakF(lmi),1))
      Ffwhm(lmi,1) = fVec(lmF_underthr(find(lmF_underthr < peakF(lmi),1,'last'))+1);
      Ffwhm(lmi,2) = fVec(end);
      Ffwhm(lmi,3) = 2*(fVec(peakF(lmi))-Ffwhm(lmi,1))+ min(diff(fVec));
    else
      Ffwhm(lmi,1) = fVec(1);
      Ffwhm(lmi,2) = fVec(end);
      Ffwhm(lmi,3) = 2*(fVec(end)-fVec(1)+min(diff(fVec)));
    end
    
    lmT_underthr = find(squeeze(TFR(peakF(lmi),:,ti) < peakpower(lmi)/2)); %Indices of TFR times of < half max power at the freq of a given local peak
    if ~isempty(find(lmT_underthr < peakT(lmi), 1)) && ~isempty(find(lmT_underthr > peakT(lmi), 1))
      Tfwhm(lmi,1) = tVec(lmT_underthr(find(lmT_underthr < peakT(lmi),1,'last'))+1);
      Tfwhm(lmi,2) = tVec(lmT_underthr(find(lmT_underthr > peakT(lmi),1,'first'))-1);
      Tfwhm(lmi,3) = Tfwhm(lmi,2)-Tfwhm(lmi,1)+ min(diff(tVec));
    elseif isempty(find(lmT_underthr < peakT(lmi),1)) && ~isempty(find(lmT_underthr > peakT(lmi),1))
      Tfwhm(lmi,1) = tVec(1);
      Tfwhm(lmi,2) = tVec(lmT_underthr(find(lmT_underthr > peakT(lmi),1,'first'))-1);
      Tfwhm(lmi,3) = 2*(Tfwhm(lmi,2)-tVec(peakT(lmi)))+ min(diff(tVec));
    elseif ~isempty(find(lmT_underthr < peakT(lmi),1)) && isempty(find(lmT_underthr > peakT(lmi),1))
      Tfwhm(lmi,1) = tVec(lmT_underthr(find(lmT_underthr < peakT(lmi),1,'last'))+1);
      Tfwhm(lmi,2) = tVec(end);
      Tfwhm(lmi,3) = 2*(tVec(peakT(lmi))-Tfwhm(lmi,1))+ min(diff(tVec));
    else
      Tfwhm(lmi,1) = tVec(1);
      Tfwhm(lmi,2) = tVec(end);
      Tfwhm(lmi,3) = 2*(tVec(end)-tVec(1)+min(diff(tVec)));
    end
  end
  
  % 12 column matrix with 1. trial index, 2. trial class, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
  % 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power, ...
  TFRlocalmax = [TFRlocalmax; ti*ones(size(peakF)) classLabels(ti)*ones(size(peakF)) fVec(peakF)' Ffwhm tVec(peakT)' Tfwhm peakpower];
  Finds_localmax = [Finds_localmax; peakF];
end

%Append 12th column: normalized local maxima peak power
medianpower = median(reshape(TFR, size(TFR,1), size(TFR,2)*size(TFR,3)), 2); %Median power at each frequency across all trials
TFRlocalmax = [TFRlocalmax TFRlocalmax(:,11)./medianpower(Finds_localmax)]; %Append column with peak power normalized to median power

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2) Pick out maxima above threshold and within the frequency band (eventBand) of interest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thr = thrFOM*medianpower; %Spectral event threshold for each frequency value
spectralEvents = TFRlocalmax((TFRlocalmax(:,3)>=eventBand(1) & TFRlocalmax(:,3)<=eventBand(2) & TFRlocalmax(:,11)>=thr(Finds_localmax)),:); %Select local maxima
clear TFRlocalmax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3) Identify and organize event features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix of event features: each row is an event
% 11 column matrix with 1. trial index, 2. hit/miss, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
% 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power
spectralEvents_columnlabel={'trialind', 'classLabels', 'maximafreq', 'lowerboundFspan', 'upperboundFspan', 'Fspan', ...
    'maximatiming', 'onsettiming', 'offsettiming', 'duration', 'maximapower', 'maximapowerFOM'};
for rci=1:numel(spectralEvents_columnlabel)
  eventsind.(spectralEvents_columnlabel{rci})=rci;
end

trialSummary.classLabels = classLabels';

trialSummary.meanpower = mean(squeeze(mean(TFR(eventBand_inds,:,:),2)) ./ repmat(medianpower(eventBand_inds),1,numTrials), 1)'; %Mean trial power normalized to frequency-specific median

suprathrTFR = TFR>=repmat(thr,1,tlength,numTrials);
trialSummary.coverage = squeeze(sum(sum(suprathrTFR(eventBand_inds,:,:),1),2)) *100 / (numel(eventBand_inds)*tlength); %Calculated in percentage

% Initialize column vectors
trialSummary.eventnumber = nan(numTrials,1);
trialSummary.meaneventpower = nan(numTrials,1);
trialSummary.meaneventduration = nan(numTrials,1);
trialSummary.meaneventFspan = nan(numTrials,1);
trialSummary.mostrecenteventtiming = nan(numTrials,1);
trialSummary.mostrecenteventpower = nan(numTrials,1);
trialSummary.mostrecenteventduration = nan(numTrials,1);
trialSummary.mostrecenteventFspan = nan(numTrials,1);

% Iterate through trials
for tri=1:numTrials
  trialSummary.eventnumber(tri)=nnz(spectralEvents(:,1)==tri);
  if nnz(spectralEvents(:,1)==tri)==0
    trialSummary.meaneventpower(tri) = 0; % traces2TFR always returns a positive value
    trialSummary.meaneventduration(tri) = 0;
    trialSummary.meaneventFspan(tri) = 0;
    trialSummary.mostrecenteventtiming(tri) = tVec(1)-mean(diff(tVec));
    trialSummary.mostrecenteventpower(tri) = 0;
    trialSummary.mostrecenteventduration(tri) = 0;
    trialSummary.mostrecenteventFspan(tri) = 0;
  else
    trialSummary.meaneventpower(tri) = mean(spectralEvents(spectralEvents(:,eventsind.trialind)==tri,eventsind.maximapowerFOM)); % traces2TFR always returns a positive value
    trialSummary.meaneventduration(tri) = mean(spectralEvents(spectralEvents(:,eventsind.trialind)==tri,eventsind.duration));
    trialSummary.meaneventFspan(tri) = mean(spectralEvents(spectralEvents(:,eventsind.trialind)==tri,eventsind.Fspan));
    trialSummary.mostrecenteventtiming(tri) = spectralEvents(find(spectralEvents(:,eventsind.trialind)==tri,1,'last'), eventsind.maximatiming);
    trialSummary.mostrecenteventpower(tri) = spectralEvents(find(spectralEvents(:,eventsind.trialind)==tri,1,'last'), eventsind.maximapowerFOM);
    trialSummary.mostrecenteventduration(tri) = spectralEvents(find(spectralEvents(:,eventsind.trialind)==tri,1,'last'), eventsind.duration);
    trialSummary.mostrecenteventFspan(tri) = spectralEvents(find(spectralEvents(:,eventsind.trialind)==tri,1,'last'), eventsind.Fspan);
  end
end

% Event dependent features (mean power, mean length, most recent timing): 
% need special treatment for zero event trials
specialFeat.field = {'meaneventpower','meaneventduration','meaneventFspan','mostrecenteventtiming',...
    'mostrecenteventpower','mostrecenteventduration','mostrecenteventFspan'};

% Percent change from mean (PCM)
trialSum_featNames = fieldnames(trialSummary);
for feat_i=2:numel(trialSum_featNames)
    pcm_name = [trialSum_featNames{feat_i},'_pcm'];
    feature = trialSummary.(trialSum_featNames{feat_i});
    
    % Control for features that need special treatment
    validtrials = trialSummary.eventnumber>0 ; %Trials that do have events
    if ismember(trialSum_featNames{feat_i},specialFeat.field)
        trialSummary.(pcm_name) = 100 * (feature-mean(feature(validtrials))) ./ repmat(abs(mean(feature(validtrials))),numTrials,1);
    else
        trialSummary.(pcm_name) = 100 * (feature-mean(feature))./repmat(abs(mean(feature)),numTrials,1);
    end
end

% Inter-event interval (IEI)
ieitemp=diff(spectralEvents(:,eventsind.maximatiming));
sametrial=(diff(spectralEvents(:,eventsind.trialind))==0);
IEI.IEI_all = ieitemp(sametrial);

for cls_i=1:numel(classes)
    fieldName = ['IEI_',num2str(classes(cls_i))];
    iei_class=diff(spectralEvents(spectralEvents(:,eventsind.classLabels)==classes(cls_i),eventsind.maximatiming));
    sametrial_class=(diff(spectralEvents(spectralEvents(:,eventsind.classLabels)==classes(cls_i),eventsind.trialind)) == 0);
    IEI.(fieldName) = iei_class(sametrial_class);
end

% Assign output structure with 3 main branches: trial-level summary 
% (TrialSummary), trial-specific events (Events), and mean inter-event 
% interval across trials (IEI)
specEv_struct.TrialSummary = struct('NumTrials',numTrials,'SpecialFeatures',specialFeat,'TrialSummary',trialSummary);
specEv_struct.Events = struct('EventBand',eventBand,'ThrFOM',thrFOM,'MedianPower',medianpower,'Threshold',thr,'Events',struct(spectralEvents_columnlabel{1},spectralEvents(:,1),...
    spectralEvents_columnlabel{2},spectralEvents(:,2),spectralEvents_columnlabel{3},spectralEvents(:,3),spectralEvents_columnlabel{4},spectralEvents(:,4),...
    spectralEvents_columnlabel{5},spectralEvents(:,5),spectralEvents_columnlabel{6},spectralEvents(:,6),spectralEvents_columnlabel{7},spectralEvents(:,7),...
    spectralEvents_columnlabel{8},spectralEvents(:,8),spectralEvents_columnlabel{9},spectralEvents(:,9),spectralEvents_columnlabel{10},spectralEvents(:,10),...
    spectralEvents_columnlabel{11},spectralEvents(:,11),spectralEvents_columnlabel{12},spectralEvents(:,12)));
specEv_struct.IEI = IEI;
end