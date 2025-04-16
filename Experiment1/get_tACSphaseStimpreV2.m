function [dataHilbert, lphpData_Stim,dataCutGaps] = get_tACSphaseStimpreV2(INdir,Session)

%1.-Defining trials
cfg              = [];
cfg.trialfun     = 'trialfun_gapDetectTACS';
cfg.headerfile   = [INdir '_tACS' num2str(Session) '.vhdr'];
cfg.datafile     = [INdir '_tACS' num2str(Session) '.eeg'];
cfg              = ft_definetrial(cfg);
trlstim_EEG      = cfg.trl;

%read the data
cfg              = [];
cfg.channel  = {'tACS1' 'tACS2' 'StimTrack'};
cfg.continuous   = 'yes';
cfg.datafile     = [INdir '_tACS' num2str(Session) '.eeg'];
data_EEG         = ft_preprocessing(cfg);

%browse through the data
cfg          = [];
cfg.channel  = {'tACS1' 'tACS2' 'StimTrack'};
cfg.viewmode = 'vertical';
ft_databrowser(cfg,data_EEG)

% 2. -bandpass filter before anything???
cfg            = [];
cfg.channel    = {'tACS1' 'tACS2'};
cfg.continuous = 'yes';
cfg.lpfilter   = 'yes';                              % apply highpass filter
cfg.lpfreq     = 10;
%cfg.hpinstabilityfix = deal with filter instability, 'no', 'reduce', 'split' (default  = 'no')
cfg.lpfiltwintype = 'hann';%highpass window type, 'hann' or 'hamming' (default) or 'blackman' or 'kaiser' (firws)
lpData_EEG     = ft_preprocessing(cfg,data_EEG);

cfg = [];
cfg.channel = 'StimTrack';
StimTrack_EEG  = ft_preprocessing(cfg,data_EEG);

cfg            = [];
cfg.channel    = {'tACS1' 'tACS2'};
cfg.continuous = 'yes';
cfg.hpfilter   = 'yes';                              % apply highpass filter
cfg.hpfreq     = 1;
cfg.hpfiltord  = 5;
lphpData_EEG     = ft_preprocessing(cfg,lpData_EEG);

if (trlstim_EEG(1,1)-1)<0 %in case the EEG started late and we cannot use 1 s before stim, in this case the trial can be excluded later
    trlstim_EEG(1,:)     = [];
    disp('......WARNING NOT ENOUGH PRE-STIM TIME FOR TRIAL 1... TRIAL will be excluded')
end

trCurrStimTrack = trlstim_EEG(trlstim_EEG(:,4)==1,:);
trCurrGapTrack = trlstim_EEG(trlstim_EEG(:,4)==2,:);

%Cut Data according to stimulus periods
cfg           = [];
cfg.trl       = trCurrStimTrack;
%cfg.trl(:,1)=cfg.trl(:,1)-cfg.trl(:,3)
lphpData_Stim = ft_redefinetrial(cfg, lphpData_EEG);
DataStimt     = ft_redefinetrial(cfg, StimTrack_EEG);

%inspect the data
cfg = [];
%cfg.channel= 'EEG';
ft_databrowser(cfg,lphpData_Stim);
ft_databrowser(cfg,DataStimt);

%Cut data according to gaps
cfg           = [];
cfg.trl       = trCurrGapTrack;
dataCutGaps   = ft_redefinetrial(cfg, lphpData_EEG);
StimTrackGaps = ft_redefinetrial(cfg, StimTrack_EEG);

%inspect the data
cfg = [];
%cfg.channel= 'EEG';
ft_databrowser(cfg,dataCutGaps);
ft_databrowser(cfg,StimTrackGaps);

%compute hilbert and get phase angle on the continuous data
dataHilbert.tACS1           = hilbert(lphpData_EEG.trial{1,1}(1,:));
dataHilbert.tACS2           = hilbert(lphpData_EEG.trial{1,1}(2,:));
dataHilbert.tACS1angle      = angle(dataHilbert.tACS1);
dataHilbert.tACS2angle      = angle(dataHilbert.tACS2);
dataHilbert.tACS1angleStim  = dataHilbert.tACS1angle((trCurrStimTrack(:,1)-trCurrStimTrack(:,3))');
dataHilbert.tACS2angleStim  = dataHilbert.tACS2angle((trCurrStimTrack(:,1)-trCurrStimTrack(:,3))');
dataHilbert.tACS1angleGap   = dataHilbert.tACS1angle((trCurrGapTrack(:,1)-trCurrGapTrack(:,3))');
dataHilbert.tACS2angleGap   = dataHilbert.tACS2angle((trCurrGapTrack(:,1)-trCurrGapTrack(:,3))');
dataHilbert.trCurrStimTrack = trCurrStimTrack;
dataHilbert.trCurrGapTrack  = trCurrGapTrack;
%plotting phase info to visualize

figure
subplot(2,2,1),circ_plot(dataHilbert.tACS1angleStim','hist')
title('tACS1 phase at Stim onset')
subplot(2,2,2),circ_plot(dataHilbert.tACS2angleStim','hist')
title('tACS2 phase at Stim onset')
subplot(2,2,3),circ_plot(dataHilbert.tACS1angleGap','hist')
title('tACS1 phase at Gap onset')
subplot(2,2,4),circ_plot(dataHilbert.tACS2angleGap','hist')
title('tACS2 phase at Gap onset')
end

