function specEv_struct = spectralevents_find(eventBand, thrFOM, tVec, fVec, TFR, classLabels)
% 1) Retrieve all local maxima in TFR using imregionalmax
% 2) Pick out maxima for the frequency band (eventBand) of interest
% 3) 


eventBand_inds = fVec(fVec>=eventBand(1) & fVec<=eventBand(2)); %Indices of freq vector within eventBand
trial_duration = tVec(end)-tVec(1); %Duration of time trial
flength = size(TFR,1); %Number of elements in discrete frequency spectrum
tlength = size(TFR,2); %Number of points in time
numtrials = size(TFR,3); %Number of trials
%time_inds = (1:length(tVec)); %Indices of timepoints corresponding to time vector

%%%% (1) %%%%

% TFRlocalmax: 12 column matrix for storing local max event metrics: trial index,
% hit/miss, maxima frequency, lowerbound frequency, upperbound frequency,
% frequency span, maxima timing, event onset timing, event offset timing, event duration, maxima power, maxima/median power
TFRlocalmax = [];

% Finds_localmax: stores peak frequecy at each local max (columns) for each trial (rows)
Finds_localmax = [];

% Retrieve all local maxima in TFR using imregionalmax
for ti=1:numtrials
  [peakF,peakT]=find(imregionalmax(squeeze(TFR(:,:,ti)))); %Indices of max local power
  peakpower=TFR(find(imregionalmax(squeeze(TFR(:,:,ti))))+(ti-1)*flength*tlength); %Power values at local maxima (vector; compiles across frequencies and time)
  %peakpower=TFR(peakF,peakT,ti); %Power values at local maxima
  
  % Find local maxima lowerbound, upperbound, and full width at half max
  % for both frequency and time
  Ffwhm=NaN(numel(peakpower),3); %2D matrix for freq-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
  Tfwhm=NaN(numel(peakpower),3); %2D matrix for time-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
  for lmi=1:numel(peakpower)
    lmF_underthr=find(squeeze(TFR(:,peakT(lmi),ti) < peakpower(lmi)/2)); %Indices of TFR frequencies of < half max power at the time of a given local peak
    if ~isempty(find(lmF_underthr < peakF(lmi), 1)) && ~isempty(find(lmF_underthr > peakF(lmi), 1))
      Ffwhm(lmi,1)=fVec(lmF_underthr(find(lmF_underthr < peakF(lmi),1,'last'))+1);
      Ffwhm(lmi,2)=fVec(lmF_underthr(find(lmF_underthr > peakF(lmi),1,'first'))-1);
      Ffwhm(lmi,3)=Ffwhm(lmi,2)-Ffwhm(lmi,1)+ min(diff(fVec));
    elseif isempty(find(lmF_underthr < peakF(lmi),1)) && ~isempty(find(lmF_underthr > peakF(lmi),1))
      Ffwhm(lmi,1)=fVec(1);
      Ffwhm(lmi,2)=fVec(lmF_underthr(find(lmF_underthr > peakF(lmi),1,'first'))-1);
      Ffwhm(lmi,3)=2*(Ffwhm(lmi,2)-fVec(peakF(lmi)))+ min(diff(fVec));
    elseif ~isempty(find(lmF_underthr < peakF(lmi),1)) && isempty(find(lmF_underthr > peakF(lmi),1))
      Ffwhm(lmi,1)=fVec(lmF_underthr(find(lmF_underthr < peakF(lmi),1,'last'))+1);
      Ffwhm(lmi,2)=fVec(end);
      Ffwhm(lmi,3)=2*(fVec(peakF(lmi))-Ffwhm(lmi,1))+ min(diff(fVec));
    else
      Ffwhm(lmi,1)=fVec(1);
      Ffwhm(lmi,2)=fVec(end);
      Ffwhm(lmi,3)=2*(fVec(end)-fVec(1)+min(diff(fVec)));
    end
    
    lmT_underthr=find(squeeze(TFR(peakF(lmi),:,ti) < peakpower(lmi)/2)); %Indices of TFR times of < half max power at the freq of a given local peak
    if ~isempty(find(lmT_underthr < peakT(lmi), 1)) && ~isempty(find(lmT_underthr > peakT(lmi), 1))
      Tfwhm(lmi,1)=tVec(lmT_underthr(find(lmT_underthr < peakT(lmi),1,'last'))+1);
      Tfwhm(lmi,2)=tVec(lmT_underthr(find(lmT_underthr > peakT(lmi),1,'first'))-1);
      Tfwhm(lmi,3)=Tfwhm(lmi,2)-Tfwhm(lmi,1)+ min(diff(tVec));
    elseif isempty(find(lmT_underthr < peakT(lmi),1)) && ~isempty(find(lmT_underthr > peakT(lmi),1))
      Tfwhm(lmi,1)=tVec(1);
      Tfwhm(lmi,2)=tVec(lmT_underthr(find(lmT_underthr > peakT(lmi),1,'first'))-1);
      Tfwhm(lmi,3)=2*(Tfwhm(lmi,2)-tVec(peakT(lmi)))+ min(diff(tVec));
    elseif ~isempty(find(lmT_underthr < peakT(lmi),1)) && isempty(find(lmT_underthr > peakT(lmi),1))
      Tfwhm(lmi,1)=tVec(lmT_underthr(find(lmT_underthr < peakT(lmi),1,'last'))+1);
      Tfwhm(lmi,2)=tVec(end);
      Tfwhm(lmi,3)=2*(tVec(peakT(lmi))-Tfwhm(lmi,1))+ min(diff(tVec));
    else
      Tfwhm(lmi,1)=tVec(1);
      Tfwhm(lmi,2)=tVec(end);
      Tfwhm(lmi,3)=2*(tVec(end)-tVec(1)+min(diff(tVec)));
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


%%%% (2) %%%%

%rhythmtfr = TFR(:,time_inds,:);
%rhythmmedian=median(reshape(rhythmtfr, size(rhythmtfr,1), size(rhythmtfr,2)*size(rhythmtfr,3)), 2);

%rhythmlmi = find(TFRlocalmax(:,3)>=fVec(eventBand_inds(1)) & TFRlocalmax(:,3)<=fVec(eventBand_inds(end)) & TFRlocalmax(:,7)>=-trial_duration & TFRlocalmax(:,7)<0);
eventBand_lmi = find(TFRlocalmax(:,3)>=eventBand(1) & TFRlocalmax(:,3)<=eventBand(2)); %Indices of local maxima within the specified event freq band

%prestimrhythmlocalmax = zeros(numel(eventBand_lmi),size(TFRlocalmax,2)+1);
%prestimrhythmlocalmax(:,size(TFRlocalmax,2)+1)=prestimrhythmlocalmax(:,11)./medianpower(Finds_localmax(eventBand_lmi)); %Add column with peak power normalized to median power
eventBand_locmax = TFRlocalmax(eventBand_lmi,:); %Select local maxima within the specified event freq band

%TFRmedian=median(reshape(TFR, size(TFR,1), size(TFR,2)*size(TFR,3)), 2);
%TFRlocalmax=[TFRlocalmax TFRlocalmax(:,11)./medianpower(TFRlocalmax(:,3))]; %Add column with peak power normalized to median power


thr = thrFOM*medianpower; %Spectral event threshold

if flength~=length(fVec) || tlength~=length(tVec) || numtrials~=length(classLabels)
    error('mismatch in input parameter dimensions')
end

%%% matrix if event features: each row is an event
%load(strcat(resultPath, rhythmid, 'localmax_', datatypeid, '_', subject_list{subjind}, '.mat'));
clear TFRlocalmax
% 11 column matrix with 1. trial index, 2. hit/miss, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
% 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power
rhythmevents_columnlabel={'trialind', 'classLabels', 'maximafreq', 'lowerboundFspan', 'upperboundFspan', 'Fspan', ...
    'maximatiming', 'onsettiming', 'offsettiming', 'duration', 'maximapower', 'maximapowerFOM'};
for rci=1:numel(rhythmevents_columnlabel)
  eventsind.(rhythmevents_columnlabel{rci})=rci;
end
eventinds = eventBand_locmax(:,eventsind.maximapower)>=thr(Finds_localmax(eventBand_lmi));
rhythmevents = zeros(nnz(eventinds),numel(rhythmevents_columnlabel));
rhythmevents(:,1:size(eventBand_locmax,2)) = eventBand_locmax(eventinds,:);

%%% trial summary of event features: each row is a trial
clear sumind
sumind.meanpower=1;
sumind.coverage=2;
sumind.eventnumber=3;
sumind.meaneventpower=4; 
sumind.meaneventduration=5;
sumind.meaneventFspan=6;
sumind.mostrecenteventtiming=7;
sumind.mostrecenteventpower=8; 
sumind.mostrecenteventduration=9;
sumind.mostrecenteventFspan=10;

indnames=fieldnames(sumind);
trialSummary=zeros(numtrials,numel(indnames));

% NORMALIZE BASED ON FREQUENCY-SPECIFIC MEDIAN
trialSummary(:,sumind.meanpower)= mean(squeeze(mean(TFR(eventBand_inds,:,:),2)) ./ repmat(medianpower(eventBand_inds),1,numtrials), 1);

suprathrTFR = TFR>=repmat(thr,1,tlength,numtrials);
trialSummary(:,sumind.coverage) = squeeze(sum(sum(suprathrTFR(eventBand_inds,:,:),1),2)) *100 / (numel(eventBand_inds)*tlength); % calculated in percentage

% iterate through trials
for tri=1:numtrials
  trialSummary(tri,sumind.eventnumber)=nnz(rhythmevents(:,1)==tri);
  if nnz(rhythmevents(:,1)==tri)==0
    trialSummary(tri,sumind.meaneventpower)=0; % traces2TFR always returns a positive value
    trialSummary(tri,sumind.meaneventduration)=0;
    trialSummary(tri,sumind.meaneventFspan)=0;
    trialSummary(tri,sumind.mostrecenteventtiming)=-trial_duration-mean(diff(tVec));
    trialSummary(tri,sumind.mostrecenteventpower)=0;
    trialSummary(tri,sumind.mostrecenteventduration)=0;
    trialSummary(tri,sumind.mostrecenteventFspan)=0;
  else
    trialSummary(tri,sumind.meaneventpower)=mean(rhythmevents(rhythmevents(:,eventsind.trialind)==tri,eventsind.maximapowerFOM)); % traces2TFR always returns a positive value
    trialSummary(tri,sumind.meaneventduration)=mean(rhythmevents(rhythmevents(:,eventsind.trialind)==tri,eventsind.duration));
    trialSummary(tri,sumind.meaneventFspan)=mean(rhythmevents(rhythmevents(:,eventsind.trialind)==tri,eventsind.Fspan));
    trialSummary(tri,sumind.mostrecenteventtiming)= rhythmevents(find(rhythmevents(:,eventsind.trialind)==tri,1,'last'), eventsind.maximatiming);
    trialSummary(tri,sumind.mostrecenteventpower)= rhythmevents(find(rhythmevents(:,eventsind.trialind)==tri,1,'last'), eventsind.maximapowerFOM);
    trialSummary(tri,sumind.mostrecenteventduration)= rhythmevents(find(rhythmevents(:,eventsind.trialind)==tri,1,'last'), eventsind.duration);
    trialSummary(tri,sumind.mostrecenteventFspan)= rhythmevents(find(rhythmevents(:,eventsind.trialind)==tri,1,'last'), eventsind.Fspan);
  end
end
    


%%% PCM and IEI

specialinds=zeros(7,1); %event depentdent parameters (mean power, mean length, most recent timing): need special treatment for zero event trials
specialinds(1)=sumind.meaneventpower;
specialinds(2)=sumind.meaneventduration;
specialinds(3)=sumind.meaneventFspan;
specialinds(4)=sumind.mostrecenteventtiming;
specialinds(5)=sumind.mostrecenteventpower;
specialinds(6)=sumind.mostrecenteventduration;
specialinds(7)=sumind.mostrecenteventFspan;
if nnz(specialinds==0)
    error('Special inds were not assigned correctly.')
end

% percent change from mean
trialSummary_pcm = 100 * (trialSummary-mean(trialSummary,1))./repmat(abs(mean(trialSummary,1)),numtrials,1);
validtrials = trialSummary(:,sumind.eventnumber)>0 ; % trials that do have events
trialSummary_pcm(:,specialinds)= 100 * (trialSummary(:,specialinds)-mean(trialSummary(validtrials,specialinds),1)) ./ repmat(abs(mean(trialSummary(validtrials,specialinds),1)),numtrials,1);
    
ieitemp=diff(rhythmevents(:,eventsind.maximatiming));
sametrial=(diff(rhythmevents(:,eventsind.trialind))==0);
rhythmIEI= ieitemp(sametrial);

ieitempY=diff(rhythmevents(rhythmevents(:,eventsind.classLabels)==1,eventsind.maximatiming));
sametrialY=(diff(rhythmevents(rhythmevents(:,eventsind.classLabels)==1,eventsind.trialind)) == 0);
rhythmIEIY= ieitempY(sametrialY);

ieitempN=diff(rhythmevents(rhythmevents(:,eventsind.classLabels)==0,eventsind.maximatiming));
sametrialN=(diff(rhythmevents(rhythmevents(:,eventsind.classLabels)==0,eventsind.trialind)) == 0);
rhythmIEIN= ieitempN(sametrialN);

% Assign output structure with 3 main branches: trial summary (TrialSummary), trial-specific events (Events), and mean inter-event interval across trials (IEI)
%specEv_struct.TFR = struct('TFR',TFR,'TVec',tVec,'FVec',fVec);
specEv_struct.TrialSummary = struct('ClassLabels',classLabels','NumTrials',numtrials,'SumInd',sumind,'SpecialInds',specialinds,'TrialSummary',trialSummary,'TrialSummary_PCM',trialSummary_pcm);
specEv_struct.EventParam = struct('EventBand',eventBand,'ThrFOM',thrFOM);
specEv_struct.Events = struct('MedianPower',medianpower,'Threshold',thr,'RhythmEvents',struct(rhythmevents_columnlabel{1},rhythmevents(:,1),...
    rhythmevents_columnlabel{2},rhythmevents(:,2),rhythmevents_columnlabel{3},rhythmevents(:,3),rhythmevents_columnlabel{4},rhythmevents(:,4),...
    rhythmevents_columnlabel{5},rhythmevents(:,5),rhythmevents_columnlabel{6},rhythmevents(:,6),rhythmevents_columnlabel{7},rhythmevents(:,7),...
    rhythmevents_columnlabel{8},rhythmevents(:,8),rhythmevents_columnlabel{9},rhythmevents(:,9),rhythmevents_columnlabel{10},rhythmevents(:,10),...
    rhythmevents_columnlabel{11},rhythmevents(:,11),rhythmevents_columnlabel{12},rhythmevents(:,12)));
specEv_struct.IEI = struct('RhythmIEI',rhythmIEI,'RhythmIEIY',rhythmIEIY,'RhythmIEIN',rhythmIEIN);
end