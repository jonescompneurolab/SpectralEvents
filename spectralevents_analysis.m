function spectralevents_analysis(specEv_struct, timeseries, TFRs, tVec, fVec)
% SPECTRALEVENT_ANALYSIS(specEv_struct,timeseries,TFRs,tVec,fVec) conducts 
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
%
% See also SPECTRALEVENTS, SPECTRALEVENTS_FIND.

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
%eventNumVec = (0:5);

figure
features = {'eventnumber','meaneventpower','meaneventduration','meaneventFspan'};
for feat_i=1:numel(features)
    feature_agg = [];
    classLabels_agg =[];
    for subj_i=1:numSubj
        featInd = specEv_struct(subj_i).TrialSummary.SumInd.(features{feat_i});
        feature_agg = [feature_agg; specEv_struct(subj_i).TrialSummary.TrialSummary(:,featInd)];
        classLabels_agg = [classLabels_agg; specEv_struct(subj_i).TrialSummary.ClassLabels];
    end
    
    classes = unique(classLabels_agg);
    [featCounts_agg,bins] = histcounts(feature_agg); %Standardize bins for all class label cases
    featProb_agg = featCounts_agg./sum(featCounts_agg);
    featProb_agg(isnan(featProb_agg)) = 0; %Correct for NaN values resulting from dividing by 0 counts
    subplot(numel(features),1,feat_i)
    
    for subj_i=1:numSubj
        %flag = false;
        featInd = specEv_struct(subj_i).TrialSummary.SumInd.(features{feat_i});
        %trials = ismember(specEv_struct(subj_i).TrialSummary.TrialSummary(:,featInd),eventNumVec); %Only consider trials with event number values within eventNumVec
        feature = specEv_struct(subj_i).TrialSummary.TrialSummary(:,featInd);
        classLabels = specEv_struct(subj_i).TrialSummary.ClassLabels;

        % Calculate probability based on event number for each type of class
        % label (rows = event number; columns = class label). Only consider
        % subjects/sessions that include trials of all class labels for each
        % event number in eventNumVec
        
        featCounts = histcounts(feature,bins);
        featProb = featCounts./sum(featCounts);
        featProb(isnan(featProb)) = 0; %Correct for NaN values resulting from dividing by 0 counts
        hold on
        plot(bins(2:end)-diff(bins)/2,featProb)
        hold off
    end
    hold on
    plot(bins(2:end)-diff(bins)/2,featProb_agg,'k-','LineWidth',2)
    hold off
        
        
        
%         for evNum_i=1:numel(eventNumVec)
%             eventNum_trials = (eventNum==eventNumVec(evNum_i));
%             feature_counts(cls_i) = histcounts(feature_PCM,'BinWidth',1);
% %             for cls_i=1:numel(classes) 
% %                 eventNum_counts(cls_i) = nnz(classLabels(eventNum_trials)==classes(cls_i)); %Counts of a given event number and given class
% %             end
% 
%             % Check to ensure all subjects/sessions include trials of a given
%             % event number across all class label types
%             if nnz(eventNum_counts==0)>0
%                 flag = true;
%                 break
%             end
% 
%             eventNum_prob = eventNum_counts./sum(eventNum_counts); %Normalize to the total counts of a given event number across classes: P(classLabel|eventNum)
%             %predLabels(eventNum_trials) = classes(find(max(eventNum_prob),1)); %Most probable class for the event number given all relevant oberservations
%         end
% 
%         % Skip subject/session if flagged
%         if flag==true
%             disp('Flag!!!!!!!')
%             AUC(subj_i) = nan;
%             continue
%         end
%         %stuff = [classLabels predLabels]
%         [x,y,~,AUC(subj_i)] = perfcurve(classLabels,eventNumPCM,'1'); %Hit probability: AUC of ROC curve
% 
%          figure
%          plot(x,y)



    %     %eventNum_prob = eventNum_counts./repmat(sum(eventNum_counts,2),1,numel(classes)); %Normalize to the total counts of a given event number across classes
    %     %eventNum_prob(isnan(eventNum_prob)) = 0; %Correct for NaN values resulting from dividing by 0
    %     predModel = double(eventNum_prob > 0.5); %Binary predictive model based on most-probable class label
    %     predLabel = nan(size(eventNum)); %Array for storing predicted class labels
    %     for trl_i=1:numTrials
    %         predLabel(trl_i) = squeeze(predModel(eventNumVec==eventNum(trl_i),2)); %Predicted class labels
    %     end
    %     
    %     prob_AUC = nan(size(eventNumVec));
    %     for evNum_i=1:numel(eventNumVec)
    %         stuff = [classLabels(eventNum==eventNumVec(evNum_i)) predLabel(eventNum==eventNumVec(evNum_i))]
    %         [~,~,~,AUC] = perfcurve(classLabels(eventNum==eventNumVec(evNum_i)),predLabel(eventNum==eventNumVec(evNum_i)),1);
    %         prob_AUC(evNum_i) = AUC; %Hit probability: AUC of ROC curve
    %     end
    % 
    %     % Calculate probability based on event power for each type of class
    %     % label
    %     
    %     
    %     %subplot(4,1,1)
    %     hold on
    %     plot(eventNumVec,prob_AUC)
    %     hold off
    %     
    %     %subplot(4,1,2)
    %     
    %    
    %     %subplot(4,1,3)
    %     
    %     
    %     %subplot(4,1,4)


end

end
