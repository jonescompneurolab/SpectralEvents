function spectralevents_analysis(specEv_struct, timeseries, TFRs, tVec, fVec)
% spectralevent_analysis(specEv_struct,timeseries,TFRs,tVec,fVec) conducts 
%   basic analysis on the spectral event features and generates spectrogram
%   and probability histogram plots comparing trial classification labels 
%   (experimental conditions or outcome states) across 
%   subjects/sessions.
%
% Inputs:
%   specEv_struct - spectralevents structure array.
%   timeseries - cell array containing time-series trials by 
%       subject/session.
%   TFRs - cell array containing time-frequency responses by
%       subject/session.
%   tVec - time vector (s) over which the time-frequency responses are 
%       shown.
%   fVec - frequency vector (Hz) over which the time-frequency responses
%       are shown.

numSubj = length(specEv_struct); %Number of subjects/sessions

% Spectrograms showing trial-by-trial events (see Figure 2 in Shin et al. eLife 2017)
for subj_i=1:numSubj
    % Extract TFR attributes for given subject/session
    TFR = TFRs{subj_i};
    classLabels = specEv_struct(subj_i).TrialSummary.ClassLabels;
    eventBand = specEv_struct(subj_i).EventParam.EventBand;
    
    % Extract event attributes for a given subject/session 
    eventThr = specEv_struct(subj_i).Events.Threshold;
    trialInd = specEv_struct(subj_i).Events.RhythmEvents.trialind;
    maximaTiming = specEv_struct(subj_i).Events.RhythmEvents.maximatiming;
    maximaFreq = specEv_struct(subj_i).Events.RhythmEvents.maximafreq;
    
    eventBand_inds = fVec(fVec>=eventBand(1) & fVec<=eventBand(2)); %Indices of freq vector within eventBand
    classes = unique(classLabels); %Array of unique class labels
    
    % Make plots for each type of class
    for cls_i=1:numel(classes)
        trial_inds = find(classLabels==classes(cls_i)); %Indices of TFR trials corresponding with the given class
        
        % Plot average TFR for a given subject/session
        avgTFR = mean(TFR(:,:,trial_inds),3);
        %figure(2*(subj_i-1)+cls_i)
        figure
        subplot(14,1,(1:4))
        %clims = [0 1200];
        imagesc([tVec(1) tVec(end)],[fVec(1) fVec(end)],avgTFR)
        %imagesc([tVec(1) tVec(end)],[fVec(1) fVec(40)],avgTFR(1:40,:),clims)
        colormap jet
        colorbar
        ylabel('Hz')
        hold on
        line(tVec',repmat(eventBand,length(tVec),1)','Color','k','LineStyle',':')
        hold off
        title(['Dataset ',num2str(subj_i),', Trial class ',num2str(classes(cls_i))])
        
        % Plot 10 randomly sampled TFR trials
        rng('default');
        randTrial_inds = randperm(numel(trial_inds),10); %Sample trial indices
        clims = [0 mean(eventThr(eventBand_inds))*1.3]; %Standardize upper spectrogram scaling limit based on the average event threshold
        for trl_i=1:10
            subplot(14,1,trl_i+4)
            imagesc([tVec(1) tVec(end)],eventBand,TFR(eventBand_inds(1):eventBand_inds(end),:,trial_inds(randTrial_inds(trl_i))),clims)
            colormap jet
            colorbar
            hold on
            plot(maximaTiming(trialInd==trial_inds(randTrial_inds(trl_i))),maximaFreq(trialInd==trial_inds(randTrial_inds(trl_i))),'w.') %Add points at event maxima
            yyaxis right
            plot(tVec,timeseries{subj_i}(:,trial_inds(randTrial_inds(trl_i))))
            hold off
        end
        xlabel('s')
    end
end

% Event feature probability histograms
eventNumTotal = [];
for subj_i=1:numSubj
    eventNumTotal = [eventNumTotal; specEv_struct(subj_i).TrialSummary.TrialSummary(:,3)];
end
eventNumVec = (0:max(eventNumTotal));

figure
for subj_i=1:numSubj
    eventNum = specEv_struct(subj_i).TrialSummary.TrialSummary(:,3);
    classLabels = specEv_struct(subj_i).TrialSummary.ClassLabels;
    classes = unique(classLabels); %Array of unique class labels
    numTrials = specEv_struct(subj_i).TrialSummary.NumTrials;

    
    % Calculate probability based on event number for each type of class
    % label (rows = event number; columns = class label)
    for cls_i=1:numel(classes)
        for evNum_i=1:numel(eventNumVec)
            eventNum_counts(evNum_i,cls_i) = nnz(classLabels(eventNum==eventNumVec(evNum_i))==classes(cls_i)); %Counts of a given event number and given class
        end
    end
    
    if numel(find(eventNum_counts==0))>0
        continue
    end
    
    eventNum_prob = eventNum_counts./repmat(sum(eventNum_counts,2),1,numel(classes)); %Normalize to the total counts of a given event number across classes
    eventNum_prob(isnan(eventNum_prob)) = 0; %Correct for NaN values resulting from dividing by 0
    predModel = double(eventNum_prob > 0.5); %Binary predictive model based on most-probable class label
    predLabel = nan(size(eventNum)); %Array for storing predicted class labels
    for trl_i=1:numTrials
        predLabel(trl_i) = squeeze(predModel(eventNum(trl_i)==eventNumVec,2)); %Predicted class labels
    end
    
    for evNum_i=1:numel(eventNumVec)
        [~,~,~,AUC] = perfcurve(classLabels(eventNum==eventNumVec(evNum_i)),predLabel(eventNum==eventNumVec(evNum_i)),1);
        %prob(evNum_i) = sum(y); %Hit probability: AUC of ROC curve
    end
    % Calculate probability based on event power for each type of class
    % label
    
    
    %subplot(4,1,1)
    hold on
    plot(eventNumVec,AUC)
    hold off
    
    %subplot(4,1,2)
    
   
    %subplot(4,1,3)
    
    
    %subplot(4,1,4)
    
    
end
end