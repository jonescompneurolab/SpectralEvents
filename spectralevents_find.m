function specEv_struct = spectralevents_find(findMethod, eventBand, thrFOM, tVec, fVec, TFR, classLabels)
% SPECTRALEVENTS_FIND Algorithm for finding and calculating spectral 
%   events on a trial-by-trial basis of a single subject/session. Uses 
%   one of three methods before further analyzing and organizing event 
%   features:
%
%   1) (Primary event detection method in Shin et al. eLife 2017): Find 
%      spectral events by first retrieving all local maxima in 
%      un-normalized TFR using imregionalmax, then selecting suprathreshold
%      peaks within the frequency band of interest. This method allows for 
%      multiple, overlapping events to occur in a given suprathreshold 
%      region and does not guarantee the presence of within-band, 
%      suprathreshold activity in any given trial will render an event.
%   2) Find spectral events by first thresholding
%      entire normalize TFR (over all frequencies), then finding local 
%      maxima. Discard those of lesser magnitude in each suprathreshold 
%      region, respectively, s.t. only the greatest local maximum in each 
%      region survives (when more than one local maxima in a region have 
%      the same greatest value, their respective event timing, freq. 
%      location, and boundaries at full-width half-max are calculated 
%      separately and averaged). This method does not allow for overlapping
%      events to occur in a given suprathreshold region and does not 
%      guarantee the presence of within-band, suprathreshold activity in 
%      any given trial will render an event.
%   3) Find spectral events by first thresholding 
%      normalized TFR in frequency band of interest, then finding local 
%      maxima. Discard those of lesser magnitude in each suprathreshold region,
%      respectively, s.t. only the greatest local maximum in each region
%      survives (when more than one local maxima in a region have the same 
%      greatest value, their respective event timing, freq. location, and 
%      boundaries at full-width half-max are calculated separately and 
%      averaged). This method does not allow for overlapping events to occur in
%      a given suprathreshold region and ensures the presence of 
%      within-band, suprathreshold activity in any given trial will render 
%      an event.
%
% specEv_struct = SPECTRALEVENTS_FIND(findMethod,eventBand,thrFOM,tVec,fVec,TFR,classLabels)
% 
% Inputs:
%   findMethod - integer value specifying which event-finding method to use
%       (1, 2, or 3). Note that the method specifies how much overlap
%       exists between events. Use 1 to replicate the method used in 
%       et al. eLife 2017.
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
%       TrialSummary (trial-level features), Events (individual event 
%       characteristics), and IEI (inter-event intervals from all trials 
%       and those associated with only a given class label).
%
% See also SPECTRALEVENTS, SPECTRALEVENTS_FIND, SPECTRALEVENTS_TS2TFR, SPECTRALEVENTS_VIS.

% Initialize general data parameters
eventBand_inds = fVec>=eventBand(1) & fVec<=eventBand(2); %Logical vector representing indices of freq vector within eventBand
if size(eventBand_inds,1)~=length(eventBand_inds)
    eventBand_inds = eventBand_inds'; %Transpose so that the dimensions correspond with the frequency-domain dimension of the TFR
end
flength = size(TFR,1); %Number of elements in discrete frequency spectrum
tlength = size(TFR,2); %Number of points in time
numTrials = size(TFR,3); %Number of trials
classes = unique(classLabels);
medianpower = median(reshape(TFR, size(TFR,1), size(TFR,2)*size(TFR,3)), 2); %Median power at each frequency across all trials
thr = thrFOM*medianpower; %Spectral event threshold for each frequency value

% Validate consistency of parameter dimensions
if flength~=length(fVec) || tlength~=length(tVec) || numTrials~=length(classLabels)
  error('Mismatch in input parameter dimensions!')
end

% Find events using the method-of-choice
spectralEvents = []; %Array for storing event results
switch findMethod
    case 1
        find_localmax_method_1;
    case 2
        find_localmax_method_2;
    case 3
        find_localmax_method_3;
    otherwise
        error('Unknown event-finding method.')
end

% Make sure this subject/session contains >1 events
if isempty(spectralEvents)
    disp('Warning!! This subject/session contains no events!!')
    specEv_struct = struct('TrialSummary',[],'Events',[],'IEI',[]);
    return;
end

% Identify and organize event features

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
trialSummary.coverage = squeeze(sum(sum(suprathrTFR(eventBand_inds,:,:),1),2)) *100 / (nnz(eventBand_inds)*tlength); %Calculated in percentage

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






    function find_localmax_method_1
    % 1st event-finding method (primary event detection method in Shin et 
    % al. eLife 2017): Find spectral events by first retrieving all local 
    % maxima in un-normalized TFR using imregionalmax, then selecting 
    % suprathreshold peaks within the frequency band of interest. This 
    % method allows for multiple, overlapping events to occur in a given 
    % suprathreshold region and does not guarantee the presence of 
    % within-band, suprathreshold activity in any given trial will render 
    % an event.

        % spectralEvents: 12 column matrix for storing local max event metrics: trial 
        % index, hit/miss, maxima frequency, lowerbound frequency, upperbound 
        % frequency, frequency span, maxima timing, event onset timing, event 
        % offset timing, event duration, maxima power, maxima/median power
        spectralEvents = [];

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
            spectralEvents = [spectralEvents; ti*ones(size(peakF)) classLabels(ti)*ones(size(peakF)) fVec(peakF)' Ffwhm tVec(peakT)' Tfwhm peakpower peakpower./medianpower(peakF)];
            Finds_localmax = [Finds_localmax; peakF];
        end

        % Pick out maxima above threshold and within the frequency band of interest
        spectralEvents = spectralEvents((spectralEvents(:,3)>=eventBand(1) & spectralEvents(:,3)<=eventBand(2) & spectralEvents(:,11)>=thr(Finds_localmax)),:); %Select local maxima
    end

    function find_localmax_method_2
    % 2nd event-finding method: Find spectral events by first thresholding
    % entire normalize TFR (over all frequencies), then finding local 
    % maxima. This method does not allow for overlapping events to occur in
    % a given suprathreshold region and does not guarantee the presence of 
    % within-band, suprathreshold activity in any given trial will render 
    % an event.

        % spectralEvents: 12 column matrix for storing local max event metrics: trial 
        % index, hit/miss, maxima frequency, lowerbound frequency, upperbound 
        % frequency, frequency span, maxima timing, event onset timing, event 
        % offset timing, event duration, maxima power, maxima/median power
        spectralEvents = [];

        % Retrieve local maxima in normalized TFR using imregionalmax, 
        % discard those of lesser (un-normalized) magnitude in each suprathreshold 
        % region, respectively, and characterize event boundaries (at half max)
        for ti=1:numTrials
            TFR_ST = squeeze(TFR(:,:,ti))./medianpower; %Suprathreshold TFR: first isolate 2D TFR matrix and normalize
            TFR_ST(TFR_ST<thrFOM) = 0; %Set infrathreshold values to zero

            % Find all local maxima in suprathreshold TFR
            TFR_LM = TFR_ST.*imregionalmax(TFR_ST); %Threshold TFR at each respective local maximum
            numTotalPeaks = nnz(TFR_LM);

            % Escape this iteration when this trial contains no suprathreshold
            % local maxima
            if numTotalPeaks==0
                continue
            end

            % Find max peak in each respective suprathreshold region
            [~,regions,numReg,~] = bwboundaries(TFR_ST>=thrFOM); %Separate suprathreshold regions
            evPeakF = cell(1,numReg);
            evPeakT = cell(1,numReg);
            evPeakpower = nan(numReg,1);
            for reg_i=1:numReg
                region = zeros(size(TFR_ST)); %Initialize a blank image that will contain a single region
                region(regions==reg_i) = 1; %Set elements (pixels) in region to the value 1
                TFR_reg = TFR_LM.*region; %Regional local maxima
                [peakF_reg,peakT_reg] = find(TFR_reg); %Indices of regional local maxima
                peakpower_reg = TFR(find(TFR_reg)+(ti-1)*flength*tlength); %Power values at regional local maxima
                maxPeakpower = max(peakpower_reg);
                maxPeak_inds = find(peakpower_reg==maxPeakpower); %Indices of all instances where local maxima have the max peak power
                evPeakF{reg_i} = peakF_reg(maxPeak_inds); %Select TFR indices at max regional peak
                evPeakT{reg_i} = peakT_reg(maxPeak_inds); %Select TFR indices at max regional peak
                evPeakpower(reg_i) = maxPeakpower(1);
            end

            % Find local maxima lowerbound, upperbound, and full width at half max
            % for both frequency and time
            evBndsF = nan(numReg,3);
            evBndsT = nan(numReg,3);
            evPeakF_inds = nan(numReg,1);
            evPeakT_inds = nan(numReg,1);
            evPeakpower_norm = nan(numReg,1);
            for reg_i=1:numReg
                numRegPeaks = numel(evPeakF{reg_i});
                peakF = evPeakF{reg_i};
                peakT = evPeakT{reg_i};
                peakpower = evPeakpower(reg_i);

                Ffwhm = nan(numRegPeaks,3); %2D matrix for freq-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
                Tfwhm = nan(numRegPeaks,3); %2D matrix for time-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
                peakpower_norm = nan(numRegPeaks,1); %Vector for storing the normalized power at each regional peak
                for lmi=1:numRegPeaks
                    lmF_underthr = find(squeeze(TFR(:,peakT(lmi),ti) < peakpower/2)); %Indices of TFR frequencies of < half max power at the time of a given local peak
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

                    lmT_underthr = find(squeeze(TFR(peakF(lmi),:,ti) < peakpower/2)); %Indices of TFR times of < half max power at the freq of a given local peak
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

                    peakpower_norm(lmi) = TFR_ST(peakF(lmi),peakT(lmi));
                end

                evBndsF(reg_i,:) = mean(Ffwhm,1);
                evBndsT(reg_i,:) = mean(Tfwhm,1);
                evPeakF_inds(reg_i) = round(mean(peakF));
                evPeakT_inds(reg_i) = round(mean(peakT));
                evPeakpower_norm(reg_i) = mean(peakpower_norm);
            end

          % 12 column matrix with 1. trial index, 2. trial class, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
          % 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power, ...
          spectralEvents = [spectralEvents; ti*ones(size(evPeakpower)) classLabels(ti)*ones(size(evPeakpower))...
              fVec(evPeakF_inds)' evBndsF tVec(evPeakT_inds)' evBndsT evPeakpower evPeakpower_norm];
        end

        % Pick out maxima within the frequency band of interest
        spectralEvents = spectralEvents((spectralEvents(:,3)>=eventBand(1) & spectralEvents(:,3)<=eventBand(2)),:); %Select local maxima
    end

    function find_localmax_method_3
    % 3rd event-finding method: Find spectral events by first thresholding 
    % normalized TFR in frequency band of interest, then finding local 
    % maxima. This method does not allow for overlapping events to occur in
    % a given suprathreshold region and ensures the presence of 
    % within-band, suprathreshold activity in any given trial will render 
    % an event.

        % spectralEvents: 12 column matrix for storing local max event metrics: trial 
        % index, hit/miss, maxima frequency, lowerbound frequency, upperbound 
        % frequency, frequency span, maxima timing, event onset timing, event 
        % offset timing, event duration, maxima power, maxima/median power
        spectralEvents = [];

        % Retrieve local maxima in normalized TFR using imregionalmax, 
        % discard those of lesser (un-normalized) magnitude in each suprathreshold 
        % region, respectively, and characterize event boundaries (at half max)
        for ti=1:numTrials
            TFR_ST = squeeze(TFR(:,:,ti))./medianpower; %Suprathreshold TFR: first isolate 2D TFR matrix and normalize
            TFR_ST(TFR_ST<thrFOM) = 0; %Set infrathreshold values to zero
            TFR_ST = TFR_ST.*eventBand_inds; %Set out-of-band values to zero

            % Find all local maxima in suprathreshold TFR
            TFR_LM = TFR_ST.*imregionalmax(TFR_ST); %Threshold TFR at each respective local maximum
            numTotalPeaks = nnz(TFR_LM);

            % Escape this iteration when this trial contains no suprathreshold
            % local maxima
            if numTotalPeaks==0
                continue
            end

            % Find max peak in each respective suprathreshold region
            [~,regions,numReg,~] = bwboundaries(TFR_ST>=thrFOM); %Separate suprathreshold regions
            evPeakF = cell(1,numReg);
            evPeakT = cell(1,numReg);
            evPeakpower = nan(numReg,1);
            for reg_i=1:numReg
                region = zeros(size(TFR_ST)); %Initialize a blank image that will contain a single region
                region(regions==reg_i) = 1; %Set elements (pixels) in region to the value 1
                TFR_reg = TFR_LM.*region; %Regional local maxima
                [peakF_reg,peakT_reg] = find(TFR_reg); %Indices of regional local maxima
                peakpower_reg = TFR(find(TFR_reg)+(ti-1)*flength*tlength); %Power values at regional local maxima
                maxPeakpower = max(peakpower_reg);
                maxPeak_inds = find(peakpower_reg==maxPeakpower); %Indices of all instances where local maxima have the max peak power
                evPeakF{reg_i} = peakF_reg(maxPeak_inds); %Select TFR indices at max regional peak
                evPeakT{reg_i} = peakT_reg(maxPeak_inds); %Select TFR indices at max regional peak
                evPeakpower(reg_i) = maxPeakpower(1);
            end

            % Find local maxima lowerbound, upperbound, and full width at half max
            % for both frequency and time
            evBndsF = nan(numReg,3);
            evBndsT = nan(numReg,3);
            evPeakF_inds = nan(numReg,1);
            evPeakT_inds = nan(numReg,1);
            evPeakpower_norm = nan(numReg,1);
            for reg_i=1:numReg
                numRegPeaks = numel(evPeakF{reg_i});
                peakF = evPeakF{reg_i};
                peakT = evPeakT{reg_i};
                peakpower = evPeakpower(reg_i);

                Ffwhm = nan(numRegPeaks,3); %2D matrix for freq-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
                Tfwhm = nan(numRegPeaks,3); %2D matrix for time-dimension event metrics with columns containing lowerbound, upperbound, and fwhm, respectively
                peakpower_norm = nan(numRegPeaks,1); %Vector for storing the normalized power at each regional peak
                for lmi=1:numRegPeaks
                    lmF_underthr = find(squeeze(TFR(:,peakT(lmi),ti) < peakpower/2)); %Indices of TFR frequencies of < half max power at the time of a given local peak
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

                    lmT_underthr = find(squeeze(TFR(peakF(lmi),:,ti) < peakpower/2)); %Indices of TFR times of < half max power at the freq of a given local peak
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

                    peakpower_norm(lmi) = TFR_ST(peakF(lmi),peakT(lmi));
                end

                evBndsF(reg_i,:) = mean(Ffwhm,1);
                evBndsT(reg_i,:) = mean(Tfwhm,1);
                evPeakF_inds(reg_i) = round(mean(peakF));
                evPeakT_inds(reg_i) = round(mean(peakT));
                evPeakpower_norm(reg_i) = mean(peakpower_norm);
            end

            % 12 column matrix with 1. trial index, 2. trial class, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
            % 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power, ...
            spectralEvents = [spectralEvents; ti*ones(size(evPeakpower)) classLabels(ti)*ones(size(evPeakpower))...
              fVec(evPeakF_inds)' evBndsF tVec(evPeakT_inds)' evBndsT evPeakpower evPeakpower_norm];
        end
    end
    
end