function spectralevents_analysis(specEv_struct, timeseries, TFRs, tVec, fVec)
% SPECTRALEVENT_ANALYSIS Conduct basic analysis for the purpose of visualizing dataset spectral event 
%   features and generates spectrogram and probability histogram plots.
%
% spectralevents_analysis(specEv_struct,timeseries,TFRs,tVec,fVec)
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
        
        % Calculate average TFR for a given subject/session and determine
        % number of trials to sample
        if numel(trial_inds)>10
            numSampTrials = 10;
            if numel(trial_inds)>1
                avgTFR = mean(TFR(:,:,trial_inds),3);
            else
                avgTFR = squeeze(TFR);
            end
        else
            numSampTrials = numel(trial_inds);
        end
        
        % Find sample trials to view
        rng('default');
        randTrial_inds = randperm(numel(trial_inds),numSampTrials); %Sample trial indices
        
        % Plot average TFR
        figure
        %pos_1 = [0.09 0.75 0.8 0.17];
        subplot('Position',[0.09 0.75 0.8 0.17])
        imagesc([tVec(1) tVec(end)],[fVec(1) fVec(end)],avgTFR)
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[fVec(1),eventBand,fVec(end)])
        ylabel('Hz')
        pos = get(gca,'position');
        colormap jet
        cb = colorbar;
        %pos([1,3]) = [0.91 0.01];
        set(cb,'position',[0.9 pos(2) 0.01 pos(4)])
        hold on
        line(tVec',repmat(eventBand,length(tVec),1)','Color','k','LineStyle',':')
        hold off
        title(['Dataset ',num2str(subj_i),', Trial class ',num2str(classes(cls_i))])
        
        % Plot 10 randomly sampled TFR trials
        clims = [0 mean(eventThr(eventBand_inds))*1.3]; %Standardize upper spectrogram scaling limit based on the average event threshold
        for trl_i=1:10
            %pos_2 = [0.09 0.75-(0.065*trl_i) 0.8 0.05];
            subplot('Position',[0.09 0.75-(0.065*trl_i) 0.8 0.05])
            imagesc([tVec(1) tVec(end)],eventBand,TFR(eventBand_inds(1):eventBand_inds(end),:,trial_inds(randTrial_inds(trl_i))),clims)
            x_ticks = get(gca,'xtick');
            x_tick_labels = get(gca,'xticklabels');
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',eventBand)
            pos = get(gca,'position');
            colormap jet
            
            % Overlay locations of event peaks and the waveform corresponding with each trial
            hold on
            plot(maximaTiming(trialInd==trial_inds(randTrial_inds(trl_i))),maximaFreq(trialInd==trial_inds(randTrial_inds(trl_i))),'w.') %Add points at event maxima
            yyaxis right
            plot(tVec,timeseries{subj_i}(:,trial_inds(randTrial_inds(trl_i))),'w')
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])
            hold off
        end
        cb = colorbar;
        %pos([1,3]) = [0.91 0.01];
        set(cb,'position',[0.9 pos(2) 0.01 pos(4)])
        set(gca,'xtick',x_ticks)
        set(gca,'xticklabel',x_tick_labels)
        xlabel('s')
    end
end

% Event feature probability histograms (see Figure 5 in Shin et al. eLife 2017)
features = {'eventnumber','maximapowerFOM','duration','Fspan'}; %Fields within specEv_struct
feature_names = {'event number','event power (FOM)','event duration (ms)','event F-span (Hz)'}; %Full names describing each field
figure
for feat_i=1:numel(features)
    feature_agg = [];
    classLabels_agg =[];
    for subj_i=1:numSubj
        % Feature-specific considerations
        if isequal(features{feat_i},'eventnumber')
            featInd = specEv_struct(subj_i).TrialSummary.SumInd.(features{feat_i});
            feature_agg = [feature_agg; specEv_struct(subj_i).TrialSummary.TrialSummary(:,featInd)];
            classLabels_agg = [classLabels_agg; specEv_struct(subj_i).TrialSummary.ClassLabels];
        else
            if isequal(features{feat_i},'duration')
                feature_agg = [feature_agg; specEv_struct(subj_i).Events.RhythmEvents.(features{feat_i}) * 1000]; %Note: convert from s->ms
            else
                feature_agg = [feature_agg; specEv_struct(subj_i).Events.RhythmEvents.(features{feat_i})];
            end
            classLabels_agg = [classLabels_agg; specEv_struct(subj_i).Events.RhythmEvents.classLabels];
        end
    end
    
    % Calculate probability of aggregate (accross subjects/sessions) and 
    % standardize bins
    [featProb_agg,bins] = histcounts(feature_agg,'Normalization','probability'); 
    
    % Correct to show left-side dropoff of histogram if applicable
    if bins(2)-(bins(2)-bins(1))/2>0
        bins = [bins(1)-(bins(2)-bins(1)),bins];
        featProb_agg = histcounts(feature_agg,bins,'Normalization','probability');
    end
    featProb_agg(isnan(featProb_agg)) = 0; %Correct for NaN values resulting from dividing by 0 counts
    
    % Calculate and plot for each subject
    subplot(numel(features),1,feat_i)
    for subj_i=1:numSubj
        % Feature-specific considerations
        if isequal(features{feat_i},'eventnumber')
            featInd = specEv_struct(subj_i).TrialSummary.SumInd.(features{feat_i});
            feature = specEv_struct(subj_i).TrialSummary.TrialSummary(:,featInd);
        else
            feature = specEv_struct(subj_i).Events.RhythmEvents.(features{feat_i});
            if isequal(features{feat_i},'duration')
                feature = feature*1000; %Convert from s->ms
            end
        end

        % Calculate probability for each subject
        featProb = histcounts(feature,bins,'Normalization','probability');
        featProb(isnan(featProb)) = 0; %Correct for NaN values resulting from dividing by 0 counts
        hold on
        plot(bins(2:end)-diff(bins)/2,featProb)
        hold off
    end
    hold on
    plot(bins(2:end)-diff(bins)/2,featProb_agg,'k-','LineWidth',2)
    hold off
    xlim([bins(2)-(bins(2)-bins(1))/2,bins(find(cumsum(featProb_agg)>=0.95,1))]) %Lower limit: smallest mid-bin; upper limit: 95% cdf interval
    xlabel(feature_names{feat_i})
    ylabel('probability')
end

end
