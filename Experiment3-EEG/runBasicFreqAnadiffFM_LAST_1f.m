function [RawMeanNorm, FMamp,phaseLockedData, fData, data1f] = runBasicFreqAnadiffFM_LAST_1f(SFdata,FMrateAll,Stim)

%PIPELINE FOR CALCULATING POWER SPECTRUM
cfg = [];

cfg.keeptrials = 'yes';
cfg.latency    = [1 19.99];%[begin end]; removebaseline and first stim sec (to exclude evoked response)
TLDATA         = ft_timelockanalysis(cfg,SFdata);
%TLDATA = SFdata;
TLDATA.trialinfo(:,3) = FMrateAll';

cfg = [];
cfg.layout      = 'acticap-64ch-standard2.mat';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, TLDATA);

%Then apply the phase locking script.
Fs = SFdata.fsample; %sampling frequency in Hz
[phaseLockedData,FreqAll] = phaselock_Trials(TLDATA,Fs);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %sanity check with the stim
% %TLDATA.time=linspace(1/44100,20*44100,882000);
% for t=1:size(TLDATA.trial,1)
%     for elec=1:size(TLDATA.trial,2)
%  %TLDATA.trial(t,elec,1:882000) = Stim(t).FreqModulator;
%      temp=resample(Stim(t).FreqModulator,Fs,44100);
% TLDATA.trial(t,elec,:) = temp(1:size(TLDATA.trial,3));
%     end
% end
% [phaseLockedData,FreqAll] = phaselock_Trials(TLDATA,Fs);
% figure, plot(squeeze(phaseLockedData(1).data.trial(:,14,:))')
% figure, plot(squeeze(phaseLockedData(2).data.trial(:,14,:))')
% figure, plot(squeeze(phaseLockedData(3).data.trial(:,14,:))')
% figure, plot(squeeze(phaseLockedData(4).data.trial(:,14,:))')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%downsample after phase-locking the trials
cfg            = [];
cfg.resamplefs = 500;  %%% define resample frequency, not all subjects had the same sampling frequency
for cond =1:length(phaseLockedData)
    phaseLockedData(cond).data  = ft_resampledata(cfg,phaseLockedData(cond).data);
end

%run old spectra analysis
TOI          = [1 18]; % cut off the first and last second, the trial included also one sec before
freqres      = 0.02;
Fs           = phaseLockedData(1).data.cfg.resamplefs;
nfft         = 90000; %get_nfft_for_res(freqres,Fs);
 F = 1; % weight of (F)ractal components of the simulated data
    O = 1; % weight of (O)scillatory components of the simulated data
     elecselec = [4 5 6 9 10 13 14 15 33]; %this are the electrodes to focuse on
   
for FMrate = 1:length(phaseLockedData)
%Do 1/f correction?

data = phaseLockedData(FMrate).data;
    
    % compute the fractal and original spectra
    cfg               = [];
    cfg.foilim        = [.5 30];
    cfg.pad           = 'nextpow2';
    cfg.tapsmofrq     = .1;
    cfg.method        = 'mtmfft';
    cfg.output        = 'fooof_aperiodic';
    data1f(FMrate).fractal = ft_freqanalysis(cfg, data);
    cfg.output        = 'pow';
    data1f(FMrate).original = ft_freqanalysis(cfg, data);

    % subtract the fractal component from the power spectrum
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2-x1';
    data1f(FMrate).oscillatory = ft_math(cfg, data1f(FMrate).fractal, data1f(FMrate).original);
    
    % original implementation by Donoghue et al. 2020 operates through the semilog-power 
    % (linear frequency, log10-power) space and transformed back into linear-linear space.
    % thus defining an alternative expression for the oscillatory component as the quotient of 
    % the power spectrum and the fractal component
    cfg               = [];
    cfg.parameter     = 'powspctrm';
    cfg.operation     = 'x2./x1';  % equivalent to 10^(log10(x2)-log10(x1))
    data1f(FMrate).oscillatory_alt = ft_math(cfg, data1f(FMrate).fractal, data1f(FMrate).original);

    % display the spectra on a log-log scale
    figure();
    subplot(1,2,1); hold on;
    plot(log(data1f(FMrate).original.freq), log(mean(data1f(FMrate).original.powspctrm(elecselec,:),1)),'k');
    plot(log(data1f(FMrate).fractal.freq), log(mean(data1f(FMrate).fractal.powspctrm(elecselec,:),1)));
    plot(log(data1f(FMrate).fractal.freq), log(mean(data1f(FMrate).oscillatory.powspctrm(elecselec,:),1)));
    xlabel('log-freq'); ylabel('log-power'); grid on;
    legend({'original','fractal','oscillatory = spectrum-fractal'},'location','southwest');
    if F~=0 && O==0
      title('pure fractal signal');
    elseif F==0 && O~=0
      title('pure oscillatory signal');
    elseif F~=0 && O~=0
      title('mixed signal');
    end
    subplot(1,2,2); hold on;
    plot(log(data1f(FMrate).original.freq), log(mean(data1f(FMrate).original.powspctrm(elecselec,:),1)),'k');
    plot(log(data1f(FMrate).fractal.freq), log(mean(data1f(FMrate).fractal.powspctrm(elecselec,:),1)));
    plot(log(data1f(FMrate).oscillatory_alt.freq), log(mean(data1f(FMrate).oscillatory_alt.powspctrm(elecselec,:),1)));
    xlabel('log-freq'); ylabel('log-power'); grid on;
    legend({'original','fractal','oscillatory = spectrum/fractal'},'location','southwest');
    title('oscillatory = spectrum / fractal');
    
    
    %get data in the format time x trial x channel
    fData(FMrate).data2Spec = shiftdim(phaseLockedData(FMrate).data.trial(:,:,phaseLockedData(FMrate).data.time >= TOI(1) & phaseLockedData(FMrate).data.time <= TOI(2)),2);
    [fData(FMrate).amp_Tr, fData(FMrate).F_Tr, fData(FMrate).xf_Tr] = spectra(fData(FMrate).data2Spec,Fs,nfft,@hann);
    fData(FMrate).Ph_Tr = angle(fData(FMrate).F_Tr);
    
    VectLengthTrAll = circ_r(fData(FMrate).Ph_Tr,[],[],2);
    
    figure('Name',['Phase lag: ' num2str(FMrate)]);
    subplot(2,2,1)
    circ_plot(fData(FMrate).Ph_Tr(145,:,14)','pretty')
    subplot(2,2,2)
    circ_plot(fData(FMrate).Ph_Tr(361,:,14)','pretty')
    subplot(2,2,3)
    circ_plot(fData(FMrate).Ph_Tr(577,:,14)','pretty')
    subplot(2,2,4)
    circ_plot(fData(FMrate).Ph_Tr(793,:,14)','pretty')
    
    figure('Name',['Vector Length FM: ' num2str(FMrate)]);
    for e=1:size(VectLengthTrAll,3)
        subplot(7,5,e),plot(fData(FMrate).xf_Tr(fData(FMrate).xf_Tr<10),VectLengthTrAll(fData(FMrate).xf_Tr<10,1,e))
        title(phaseLockedData(FMrate).data.label{e})
    end
    subplot(7,5,34),plot(fData(FMrate).xf_Tr(fData(FMrate).xf_Tr<10),mean(VectLengthTrAll(fData(FMrate).xf_Tr<10,1,elecselec),3))
    title('Average across elecselec')
    
    subplot(7,5,35),plot(fData(FMrate).xf_Tr(fData(FMrate).xf_Tr<10),mean(VectLengthTrAll(fData(FMrate).xf_Tr<10,1,:),3))
    title('Average across ALL electrodes')
    
    phaseLockedData(FMrate).data.avg = mean(phaseLockedData(FMrate).data.trial);
    phaseLockedData(FMrate).data.avg = shiftdim(phaseLockedData(FMrate).data.avg,2);
    [fData(FMrate).amp, fData(FMrate).F, fData(FMrate).xf] = spectra(phaseLockedData(FMrate).data.avg(phaseLockedData(FMrate).data.time >= TOI(1) & phaseLockedData(FMrate).data.time <= TOI(2),1,:),Fs,nfft,@hann);
    
    figure('Name',['Amp FM: ' num2str(FMrate)]);
    for e=1:size(VectLengthTrAll,3)
        subplot(7,5,e),
        plot(fData(FMrate).xf(fData(FMrate).xf<10),mean(fData(FMrate).amp(fData(FMrate).xf<10,1,e),2))
        title(phaseLockedData(FMrate).data.label{e})
        ylabel('Amp')
        xlabel('Frequency')
        %    ylim([0 0.8])
    end
    subplot(7,5,34),plot(fData(FMrate).xf(fData(FMrate).xf<10),mean(fData(FMrate).amp(fData(FMrate).xf<10,1,elecselec),3))
    title('Average across elecselec')
    
    subplot(7,5,35),plot(fData(FMrate).xf(fData(FMrate).xf<10),mean(fData(FMrate).amp(fData(FMrate).xf<10,1,:),3))
    title('Average across ALL electrodes')
    
    frq1 = round(fData(FMrate).xf,4);
    frq2 = round(FreqAll,4);
    FMamp(FMrate,:) = mean(fData(FMrate).amp(fData(FMrate).xf<=10,1,elecselec),3)';
    freqFFT = round(fData(FMrate).xf(fData(FMrate).xf<=10),2);
    fData(FMrate).VL = VectLengthTrAll;
    FMvect(FMrate,:) = mean(fData(FMrate).VL(fData(FMrate).xf<=10,1,elecselec),3)';
    
end

%PLOT amplit freq normalized by amp at other freq
for freqofI = 1:length(FreqAll)
    [minH,freqHarm1] = min(abs(freqFFT - FreqAll(freqofI)*2));
    [minM,freq1rst]  = min(abs(freqFFT - FreqAll(freqofI)));
    disp(freqFFT(freqHarm1))
    disp(freqFFT(freq1rst))
    
    RawMeanNorm(freqofI,1:2) = [FMamp(freqofI,freq1rst) mean([FMamp(1:freqofI-1,freq1rst);FMamp(freqofI+1:end,freq1rst)])];
    
    amp1rst = FMamp(:,freq1rst);
    ampHarm1 = FMamp(:,freqHarm1);
    amp1rst(freqofI) = [];
    ampHarm1(freqofI) = [];
    RawMeanNormH(freqofI,1:2) = [mean([FMamp(freqofI,freq1rst),FMamp(freqofI,freqHarm1)]) mean([amp1rst;ampHarm1])];
    
end
RawMeanNorm(:,3) = RawMeanNorm(:,1)./RawMeanNorm(:,2);
RawMeanNormH(:,3) = RawMeanNormH(:,1)./RawMeanNormH(:,2);

figure
subplot(1,2,1)
plot(frq2,RawMeanNorm(:,3))
hold on
plot(frq2,RawMeanNorm(:,1))
plot(frq2,RawMeanNorm(:,2))
xlabel('Frequency')
ylabel('Amplitude')
legend({'normalized' 'Raw' 'Baseline'})

subplot(1,2,2)
plot(frq2,RawMeanNormH(:,3))
hold on
plot(frq2,RawMeanNormH(:,1))
plot(frq2,RawMeanNormH(:,2))
xlabel('Frequency')
ylabel('Amplitude Freq1+Harm1')
legend({'normalized' 'Raw' 'Baseline'})

%SAME FOR VECTOR LENGTH
for freqofI = 1:length(FreqAll)
    [minH,freqHarm1] = min(abs(freqFFT - FreqAll(freqofI)*2));
    [minM,freq1rst]  = min(abs(freqFFT - FreqAll(freqofI)));
    disp(freqFFT(freqHarm1))
    disp(freqFFT(freq1rst))
    
    RawMeanNormV(freqofI,1:2) = [FMvect(freqofI,freq1rst) mean([FMvect(1:freqofI-1,freq1rst);FMvect(freqofI+1:end,freq1rst)])];
    
    V1rst = FMvect(:,freq1rst);
    VHarm1 = FMvect(:,freqHarm1);
    V1rst(freqofI) = [];
    VHarm1(freqofI) = [];
    RawMeanNormHV(freqofI,1:2) = [mean([FMvect(freqofI,freq1rst),FMvect(freqofI,freqHarm1)]) mean([V1rst;VHarm1])];
    
end
RawMeanNormV(:,3) = RawMeanNormV(:,1)./RawMeanNormV(:,2);
RawMeanNormHV(:,3) = RawMeanNormHV(:,1)./RawMeanNormHV(:,2);
fData(1).RawMeanNormV = RawMeanNormV;
fData(1).RawMeanNorm = RawMeanNorm;
fData(1).RawMeanNormHV = RawMeanNormHV;
fData(1).RawMeanNormH = RawMeanNormH;

figure
subplot(1,2,1)
plot(frq2,RawMeanNormV(:,3))
hold on
plot(frq2,RawMeanNormV(:,1))
plot(frq2,RawMeanNormV(:,2))
xlabel('Frequency')
ylabel('VL')
legend({'normalized' 'Raw' 'Baseline'})

subplot(1,2,2)
plot(frq2,RawMeanNormHV(:,3))
hold on
plot(frq2,RawMeanNormHV(:,1))
plot(frq2,RawMeanNormHV(:,2))
xlabel('Frequency')
ylabel('VL Freq1+Harm1')
legend({'normalized' 'Raw' 'Baseline'})
end


