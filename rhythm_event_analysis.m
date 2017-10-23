function rhythm_event_analysis(rhythmid, datatypeid, resultPath, subject_list, subjind, thrFOM, prestim_tVec, fVec, prestim_TFR_yes_no, YorN)

thr=thrFOM*rhythmmedian;
flength=size(prestim_TFR_yes_no,1);
tlength=size(prestim_TFR_yes_no,2);
numtrials=size(prestim_TFR_yes_no,3);

if flength~=length(fVec) || tlength~=length(prestim_tVec) || numtrials~=length(YorN)
    error('mismatch in input parameter dimensions')
end

%% matrix if event features: each row is an event
load(strcat(resultPath, rhythmid, 'localmax_', datatypeid, '_', subject_list{subjind}, '.mat'));
clear TFRlocalmax
% 11 column matrix with 1. trial index, 2. hit/miss, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
% 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power
rhythmevents_columnlabel={'trialind', 'yn', 'maximafreq', 'lowerboundFspan', 'upperboundFspan', 'Fspan', ...
    'maximatiming', 'onsettiming', 'offsettiming', 'duration', 'maximapower', 'maximapowerFOM'};
for rci=1:numel(rhythmevents_columnlabel)
  eventsind.(rhythmevents_columnlabel{rci})=rci;
end
eventinds = prestimrhythmlocalmax(:,eventsind.maximapower)>=thr(Finds_localmax(rhythmlmi));
rhythmevents=zeros(nnz(eventinds),numel(rhythmevents_columnlabel));
rhythmevents(:,1:size(prestimrhythmlocalmax,2))=prestimrhythmlocalmax(eventinds,:);

%% trial summary of event features: each row is a trial
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
    trialSummary(:,sumind.meanpower)= mean(squeeze(mean(prestim_TFR_yes_no(rhythm_band_inds,:,:),2)) ./ repmat(rhythmmedian(rhythm_band_inds),1,numtrials), 1);

  suprathrTFR = prestim_TFR_yes_no>=repmat(thr, 1,tlength,numtrials);
trialSummary(:,sumind.coverage) = squeeze(sum(sum(suprathrTFR(rhythm_band_inds,:,:),1),2)) *100 / (numel(rhythm_band_inds)*tlength); % calculated in percentage

% iterate through trials
for tri=1:numtrials
  trialSummary(tri,sumind.eventnumber)=nnz(rhythmevents(:,1)==tri);
  if nnz(rhythmevents(:,1)==tri)==0
    trialSummary(tri,sumind.meaneventpower)=0; % traces2TFR always returns a positive value
    trialSummary(tri,sumind.meaneventduration)=0;
    trialSummary(tri,sumind.meaneventFspan)=0;
    trialSummary(tri,sumind.mostrecenteventtiming)=-prestimduration-mean(diff(prestim_tVec));
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
    


%% PCM and IEI

specialinds=zeros(7,1); %event depentdent parameters (mean power, mean length, most recent timing): need special treatment for zero event trials
specialinds(1)=sumind.meaneventpower;
specialinds(2)=sumind.meaneventduration;
specialinds(3)=sumind.meaneventFspan;
specialinds(4)=sumind.mostrecenteventtiming;
specialinds(5)=sumind.mostrecenteventpower;
specialinds(6)=sumind.mostrecenteventduration;
specialinds(7)=sumind.mostrecenteventFspan;
if nnz(specialinds==0)
    error('special inds were not assigned correctly.')
end

% percent change from mean
trialSummary_pcm = 100 * (trialSummary-mean(trialSummary,1))./repmat(abs(mean(trialSummary,1)),numtrials,1);
validtrials = trialSummary(:,sumind.eventnumber)>0 ; % trials that do have events
trialSummary_pcm(:,specialinds)= 100 * (trialSummary(:,specialinds)-mean(trialSummary(validtrials,specialinds),1)) ./ repmat(abs(mean(trialSummary(validtrials,specialinds),1)),numtrials,1);
    
ieitemp=diff(rhythmevents(:,eventsind.maximatiming));
sametrial=(diff(rhythmevents(:,eventsind.trialind))==0);
rhythmIEI= ieitemp(sametrial);

ieitempY=diff(rhythmevents(rhythmevents(:,eventsind.yn)==1,eventsind.maximatiming));
sametrialY=(diff(rhythmevents(rhythmevents(:,eventsind.yn)==1,eventsind.trialind)) == 0);
rhythmIEIY= ieitempY(sametrialY);

ieitempN=diff(rhythmevents(rhythmevents(:,eventsind.yn)==0,eventsind.maximatiming));
sametrialN=(diff(rhythmevents(rhythmevents(:,eventsind.yn)==0,eventsind.trialind)) == 0);
rhythmIEIN= ieitempN(sametrialN);
    

%% save stuff
save(strcat(resultPath, rhythmid, 'TrialSummary_', datatypeid, '_', subject_list{subjind}, '.mat'), 'YorN', 'numtrials', 'sumind', 'specialinds', 'trialSummary', 'trialSummary_pcm');
save(strcat(resultPath, rhythmid, 'events_', datatypeid, '_', subject_list{subjind}, '.mat'), 'rhythmmedian', 'thr', 'rhythmevents_columnlabel', 'rhythmevents');
save(strcat(resultPath, rhythmid, 'IEI_', datatypeid, '_', subject_list{subjind}, '.mat'), 'rhythmmedian', 'thr', 'rhythmIEI', 'rhythmIEIY', 'rhythmIEIN');


end