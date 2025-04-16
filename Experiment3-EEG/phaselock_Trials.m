function [phaseLockedData, FreqAll] = phaselock_Trials(avg,Fs)
avgAll       = avg;
FreqAll      = unique(avgAll.trialinfo(:,3));
numphasebins = 8;
phaseLockedData = struct([]);

for condFreq = 1:length(FreqAll) %FIRST LOOP THROUGH FREQ CONDITIONS
   
    freq       = FreqAll(condFreq);
    cfg        = [];
    cfg.trials = avgAll.trialinfo(:,3)==freq;
    avg        = ft_selectdata(cfg,avgAll);

    numsamples2shift = round(((Hz2ms(freq)/1000)*Fs)/numphasebins);

    temp.time       = avg.time(1:length(avg.time)-numsamples2shift*(numphasebins));
    temp.avg        = zeros(length(avg.label),length(temp.time));
    temp.trial      = zeros(size(avg.trial,1),size(avg.trial,2),length(temp.time));
    temp.var        = zeros(size(temp.avg));
    temp.dof        = zeros(size(temp.avg));
    temp.label      = avg.label;
    temp.trialinfo  = avg.trialinfo;
  
    for k = 1:size(avg.trial,1)
        phase = avg.trialinfo(k,2);
        for p = 1:numphasebins
            if phase==p
                temp.trial(k,:,:) = avg.trial(k,:,1+numsamples2shift*(numphasebins-p+1):length(temp.time)+numsamples2shift*(numphasebins-p+1));
            end
        end
    end
    temp.avg = squeeze(mean(temp.trial,1,'omitnan'));
    temp.var = squeeze(var(temp.trial,1,'omitnan'));
    temp.dof = repmat(size(temp.trial,1),size(temp.avg));
    phaseLockedData(condFreq).data = temp;
end
end
