function rhythm_localmax_analysis(resultPath, datatypeid, subject_list, subjind, rhythmid, rhythm_band_inds, tVec, fVec, prestim_inds, TFR_yes_no, YorN)

% TFRlocalmax: first, retrieve all local maxima in TFR using imregionalmax
% prestimrhythmlocalmax: then, pick out maxima for the frequency band of interest in the prestimulus period

prestimduration=1000; %ms
flength=size(TFR_yes_no,1);
tlength=size(TFR_yes_no,2);
numtrials=size(TFR_yes_no,3);

TFRlocalmax=[];
Finds_localmax=[];

for ti=1:numtrials
  [peakF,peakT]=find(imregionalmax(squeeze(TFR_yes_no(:,:,ti))));
  %     if ~isempty(peakT)
  peakpower=TFR_yes_no(find(imregionalmax(squeeze(TFR_yes_no(:,:,ti))))+(ti-1)*flength*tlength);
  
  Ffwhm=NaN(numel(peakpower),3); %in order of lowerbound, upperbound and fwhm
  Tfwhm=NaN(numel(peakpower),3); %in order of lowerbound, upperbound and fwhm
  % calculate local maxima frequency full width at half max
  for lmi=1:numel(peakpower)
    lmF_underthr=find(squeeze(TFR_yes_no(:,peakT(lmi),ti)<peakpower(lmi)/2));
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
    
    lmT_underthr=find(squeeze(TFR_yes_no(peakF(lmi),:,ti)<peakpower(lmi)/2));
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
  
  % % % 12 column matrix with 1. trial index, 2. hit/miss, 3. maxima frequency, 4. lowerbound frequency, 5. upperbound frequency, 6. frequency span, ...
  % % % 7. maxima timing, 8. event onset timing, 9. event offset timing, 10. event duration, 11. maxima power, 12. maxima/median power, ...
  % % {'trial index', 'yes/no', 'maxima frequency', 'lowerbound F-span', 'upperbound F-span', 'F-span', ...
  % %     'maxima timing', 'onset timing', 'offset timing', 'duration', 'maxima power', 'maxima/median power'};
  TFRlocalmax=[TFRlocalmax; ti*ones(size(peakF)) YorN(ti)*ones(size(peakF)) fVec(peakF)' Ffwhm tVec(peakT)' Tfwhm peakpower];
  Finds_localmax=[Finds_localmax; peakF];
end


rhythmtfr=TFR_yes_no(:,prestim_inds,:);
rhythmmedian=median(reshape(rhythmtfr, size(rhythmtfr,1), size(rhythmtfr,2)*size(rhythmtfr,3)), 2);
rhythmlmi = find(TFRlocalmax(:,3)>=fVec(rhythm_band_inds(1)) & TFRlocalmax(:,3)<=fVec(rhythm_band_inds(end)) & TFRlocalmax(:,7)>=-prestimduration & TFRlocalmax(:,7)<0);
prestimrhythmlocalmax = zeros(numel(rhythmlmi),size(TFRlocalmax,2)+1);
prestimrhythmlocalmax(:,1:size(TFRlocalmax,2)) = TFRlocalmax(rhythmlmi,:);

prestimrhythmlocalmax(:,size(TFRlocalmax,2)+1)=prestimrhythmlocalmax(:,11)./rhythmmedian(Finds_localmax(rhythmlmi));

TFRmedian=median(reshape(TFR_yes_no, size(TFR_yes_no,1), size(TFR_yes_no,2)*size(TFR_yes_no,3)), 2);
TFRlocalmax=[TFRlocalmax TFRlocalmax(:,11)./rhythmmedian(TFRlocalmax(:,3))];

save(strcat(resultPath, rhythmid, 'localmax_', datatypeid, '_', subject_list{subjind}, '.mat'), 'rhythm_band_inds', 'prestimduration', 'rhythmmedian', 'TFRlocalmax', 'prestimrhythmlocalmax', 'Finds_localmax', 'rhythmlmi');


end