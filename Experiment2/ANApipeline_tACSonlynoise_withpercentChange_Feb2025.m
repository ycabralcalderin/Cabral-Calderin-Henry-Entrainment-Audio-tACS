%pipeline for analysing gap detection in noise-stimulus+tACS data
%Cabral-Calderin Nov. 2022

clear
close all
clc

%DEFINE DIRECTORIES
MAINdir  = uigetdir (); %main path to raw data, e.g., '/Users/yuranny.cabral/Documents/gapDetectionProject/tACS/DATA'
BEHdir   = fullfile(MAINdir, 'Behavioral'); %raw beh data
EEGdir   = fullfile(MAINdir, 'tACS'); %raw eeg data (for tACS signal)
OUTdir   = fullfile(MAINdir, 'ANA'); %main output directory
addpath('/mnt/beegfs/users/yuranny.cabral/2019-0226-relentrain/MATLAB_toolboxes/CircStat2012a')
addpath('/mnt/beegfs/users/yuranny.cabral/2019-0226-relentrain/MATLAB_toolboxes/fieldtrip-20201128')

addpath(genpath('/mnt/beegfs/users/yuranny.cabral/2019-0226-relentrain/AnalysisScripts'))

%define subject and protocol info
SUBJlist = {'ABI24' 'ATE26' 'BAO30' 'BBB19' 'BWC24' 'CSA27' 'EDI31' 'EGA23' 'EOE26' 'ESA20' 'HAE30' 'KWA10' 'LCS16' 'LPA03' 'MDS07' 'MUL19' 'OLDO03' 'PLU28' 'SFE21' 'SKN24' 'SRL03' 'TDN16' 'VEI09' 'ZJE04'};
%EEGprefereFreList
%%%% IN EEG PREFERED FREQ EXP:
%ABI24=Bei24
%EDI31=Ei31
%LPA03=LA03
%LCS16=lys16
%SRL03=SUL03
%BWC24=WEC24

%subjects = {'ATE26' 'BAO30' 'Bei24' 'CSA27' 'EGA23' 'Ei31' 'ESA20' 'LA03' 'lys16' 'MDS07' 'SUL03' 'WEC24'};

%NOTES
% SUBJlist = {'ECY06'};
%define some general variables
NrBlocks     = 5; %in general, each Session had 5 blocks
NrTrlinBlock = 32; % 32 trials were presented within the block
NrTrAll      = NrBlocks*NrTrlinBlock; % total number of trials
NrGapBins    = 6; %number of tACS phase bins at gap onset
nBins        = NrGapBins;
ft_defaults

%this part works with the raw signal, does some preprocessing and outputs
%all the info needed for further analysis, ie.e., tACS phase and detection
%response for each gap and tACS condition.
for subj= 1:length(SUBJlist)%loop across subjects 
    % close all
    currSUBJ = SUBJlist{subj}; % current subject
    fprintf('%s%s/n','... working on subject ', currSUBJ)
    %get all files for the current participant
    files = dir([fullfile(EEGdir,currSUBJ) '_*.eeg']);
    
    %load stim file to check that the modulation depth is correct
    load(fullfile(BEHdir,[currSUBJ '/3/MAIN/STIM'],[currSUBJ '_Stim_Main_orig.mat']));
    if isfield(Stim(1).stiminfo,'ModDepth')
        if Stim(1).stiminfo.ModDepth~=0
            error 'WARNING: Modulation depth is not 0'
        end
    else
    end
    
    curroDIR = fullfile(OUTdir,currSUBJ); %subject dir for saving data 
    %WORK on the .eeg files to get the tACS signal and triggers and define
    %tACS at stimulus time
    [dataHilbert, lphpData_Stim,dataCutGaps] = get_tACSphaseStimpre_noise(fullfile(EEGdir,[currSUBJ '_tACS3']));
    
    %first get max amplitude to decide on the threshold for detecting the tACS artifact
    MaxAmp = nan(size(dataCutGaps.trial));
    MinAmp = MaxAmp;
    range  = MaxAmp;
    
    for tr=1:length(dataCutGaps.trial)
        MaxAmp(tr) = max(dataCutGaps.trial{1,tr}(1,:));
        MinAmp(tr) = min(dataCutGaps.trial{1,tr}(1,:));
        range(tr)  = MaxAmp(tr)-MinAmp(tr);
    end
    ThrAmp = min(range)+((median(range)-min(range))./5);
    
    BlockLim = 1:136:length(dataCutGaps.trial);

    blockmeanrange = nan(size(BlockLim));
    if length(dataCutGaps.trial)-BlockLim(end)<100%in case there was some problem and there were some extra trials
        BlockLim(end)=[];
    end
    for bl = 1:length(BlockLim)-1
        blockmeanrange(bl) = mean(range(BlockLim(bl):BlockLim(bl+1)-1));
    end
    blockmeanrange(bl+1) = mean(range(BlockLim(end):length(dataCutGaps.trial)));
    [minrange,idx] = min(blockmeanrange);%get block with smallest range as the sham
    
    %firstSham = str2double(inputdlg('first sham trial')); %select the first sham trial based on the channel data
    firstSham = BlockLim(idx);
    disp(firstSham)

    ThrAmpShn  = 30;%min(range)+((median(range)-min(range))./5)
    NewStimTr  = range>=ThrAmp;%zeros(1,length(dataCutGaps.trial));
    NewShamTrt = find(range<ThrAmpShn);%zeros(1,length(dataCutGaps.trial));
    NewShamTrt = NewShamTrt(NewShamTrt>firstSham & NewShamTrt<(firstSham+175)); %take only the trials within the sham block
    NewShamTr  = zeros(size(NewStimTr));
    NewShamTr(NewShamTrt)=1;
    NewShamTr = NewShamTr==1;
    if sum(NewShamTr)<122 || sum(NewShamTr)>175
        disp(sum(NewShamTr))
        error 'Sham trials were not correctly identified'
    end
    %do some plotting for visualizing
    
    figure, scatter(1:length(range),range)
    hold on
    title(currSUBJ)
    plot([1 length(range)],[ThrAmp ThrAmp])
    plot([1 length(range)],[ThrAmpShn ThrAmpShn])
    scatter(find(NewShamTr==1),range(NewShamTr==1))
    scatter(find(NewStimTr==1),range(NewStimTr==1))
    
    save (fullfile(curroDIR,'preproctACS_onlyNoise_orig'),'dataCutGaps','dataHilbert','lphpData_Stim','NewShamTr','NewStimTr','-v7.3') 
end

%%%WORKING ON ALREADY PREPROCESSED DATA
dataGroup   = struct([]);
model       = struct([]);
RTerror     = nan(1,length(length(SUBJlist)));
GroupMeanHR = nan(length(length(SUBJlist)),5);
GroupMeanRT = nan(length(length(SUBJlist)),5);
        
%parameters for the cosine FIT
X         = 0:2*pi/nBins:2*pi-2*pi/nBins; %this is my phase vector
X         = X+2*pi/nBins/2;
XTime     = 0:2*pi/nBins:2*pi-2*pi/nBins; %this is my phase vector
XTime     = XTime+2*pi/nBins/2;
lag       = pi;
intercept = 0.5;
amp       = 0.5;
params    = [lag intercept amp];

for subj= 1:length(SUBJlist)%loop across subjects %14 19

    % close all
    currSUBJ = SUBJlist{subj}; % current subject
    fprintf('%s%s/n','... working on subject ', currSUBJ)
    curroDIR = fullfile(OUTdir,currSUBJ); %subject dir for saving data

    load (fullfile(curroDIR,'preproctACS_onlyNoise_orig'),'dataCutGaps','dataHilbert','lphpData_Stim','NewShamTr','NewStimTr') 

    dataHilbert.trCurrGapTrack(:,12) = dataHilbert.trCurrGapTrack(:,6); %save RT info in column 12 to not mess up with the definitions of the otehr variables
    
    %GET TIME WITHIN SOUND and tACS Freq FOR EACH GAP
    %REMOVE GAPS THAT START BEFORE THE FIRST RECORDED STIM (for a couple of participants the recording started late)
    extraGaps = dataHilbert.trCurrGapTrack(:,1)< dataHilbert.trCurrStimTrack(1,1);
    fprintf('%s%s%s',num2str(sum(extraGaps)), ' gap(s) before first stim in Subject  ', currSUBJ)
    dataHilbert.trCurrGapTrack(extraGaps,:) = [];
    NewShamTr(extraGaps) = [];
    NewStimTr(extraGaps) = [];
    dataHilbert.tACS1angleGap (extraGaps) = [];
    dataHilbert.tACS2angleGap (extraGaps) = [];  
    
    if strcmp(currSUBJ,'LPA03') %this subject had one gap without starting stim in the middle of the session
        extraGaps = 276;
        dataHilbert.trCurrGapTrack(extraGaps,:) = [];
        NewShamTr(extraGaps) = [];
        NewStimTr(extraGaps) = [];
        dataHilbert.tACS1angleGap (extraGaps) = [];
        dataHilbert.tACS2angleGap (extraGaps) = [];
        dataHilbert.freqBlock(161) = dataHilbert.freqBlock(160);
    elseif strcmp(currSUBJ,'SFE21')
        dataHilbert.freqBlock(161) = dataHilbert.freqBlock(160);
    end
    
    for g =1:length(dataHilbert.trCurrGapTrack)
        cStim = find(dataHilbert.trCurrStimTrack(:,1) < dataHilbert.trCurrGapTrack(g,1) & dataHilbert.trCurrStimTrack(:,2) > dataHilbert.trCurrGapTrack(g,2));%find stim time from the current trial
        dataHilbert.trCurrGapTrack(g,6)  = ((dataHilbert.trCurrGapTrack(g,1)-(dataHilbert.trCurrGapTrack(g,3)))-(dataHilbert.trCurrStimTrack(cStim,1)-dataHilbert.trCurrStimTrack(cStim,3)))/dataCutGaps.fsample; %get gap onset within the trial in seconds
        dataHilbert.trCurrGapTrack(g,7)  = dataHilbert.freqBlock(cStim);%tACS frequency for the trial
        dataHilbert.trCurrGapTrack(g,10) = (dataHilbert.trCurrGapTrack(g,1)-(dataHilbert.trCurrGapTrack(g,3)))-(dataHilbert.trCurrStimTrack(cStim,1)-dataHilbert.trCurrStimTrack(cStim,3));%gap onset in samples
    end
    dataHilbert.trCurrGapTrack(:,8) = nan;
    dataHilbert.trCurrGapTrack(NewStimTr,8) = dataHilbert.trCurrGapTrack(NewStimTr,7); %add frequency information only for the trials with either full stim or full sham
    dataHilbert.trCurrGapTrack(NewShamTr,8) = dataHilbert.trCurrGapTrack(NewShamTr,7);
        
    nBins   = NrGapBins;
    freqAll = unique(dataHilbert.trCurrGapTrack(:,7));
    
    for f = 1:length(freqAll) %loop through frequency conditions
        currGaps = dataHilbert.trCurrGapTrack(dataHilbert.trCurrGapTrack(:,8) == freqAll(f),:);
        currGaps(:,9) = dataHilbert.tACS1angleGap(dataHilbert.trCurrGapTrack(:,8) == freqAll(f))';
        GroupMeanHR(subj,f) = mean(currGaps(:,5));
        GroupMeanRT(subj,f) = nanmean(currGaps(:,12));
        figure(13)
        subplot(1,5,f)
        circ_plot(currGaps(:,9),'hist')
        
        %Group by Hypothetical oscillator at different frequencies to
        %simulate rhyhtmic sampling
        
        fs = lphpData_Stim.fsample; % Sampling rate 1000 Hz = 1 kHz
        t  = 0:1/fs:max(currGaps(:,6))+1; % take a sample of about 20 seconds
        binEdges = linspace(-pi,pi,nBins+1);
        dataGroup(subj).currGaps = currGaps;
        freqImg = 0.1:.1:8;
        for c = 1:length(freqImg)
            omega  = 2*pi*freqImg(c); % omega 1
            Stim2  = cos(omega.*t);
            phase2 = wrapToPi(omega.*t);
            currGaps(:,11) = phase2(currGaps(:,10)); %phase at hypothetical cos
            dataGroup(subj).currGaps(:,13+c) = currGaps(:,11);%saving phases for the imaginary cos
            [dataGroup(subj).countImgCos(f,:,c),EDGES1,binIdx] = histcounts(currGaps(:,11),binEdges);          
            dataGroup(subj).HR_ImgCosHit(f,:,c) = zeros(1,nBins,1);
            dataGroup(subj).HR_ImgCosRT(f,:,c) = zeros(1,nBins,1);
          
            for bin = 1:nBins
                if any(binIdx==bin)
                  %  figure, circ_plot(currGaps(binIdx==bin,11))
                    dataGroup(subj).HR_ImgCosHit(f,bin,c)   = mean(currGaps(binIdx==bin,5)); %Hit Rate
                    dataGroup(subj).HR_ImgCosRT(f,bin,c)    = nanmean(currGaps(binIdx==bin,12)); %RT
                end
            end
        end
        if f>1 %no tACS phase can be estimated from sham condition because in this case sham=no stim
        %Group by tACS phase
        [dataGroup(subj).countPh(f,:),EDGES2,binIdx] = histcounts(currGaps(:,9),binEdges);  
        dataGroup(subj).HR_PhaseBinHit(f,:) = zeros(1,nBins);
        dataGroup(subj).HR_PhaseBinRT(f,:)  = zeros(1,nBins);

        for bin = 1:nBins
            if any(binIdx==bin)
                dataGroup(subj).HR_PhaseBinHit(f,bin) = mean(currGaps(binIdx==bin,5));
                dataGroup(subj).HR_PhaseBinRT(f,bin)  = nanmean(currGaps(binIdx==bin,12));
%                 figure,circ_plot(currGaps(binIdx==bin,9))
            end
        end
       
        %fit by tACS phase Hit
        data2fitHit = [dataGroup(subj).HR_PhaseBinHit(f,:) dataGroup(subj).HR_PhaseBinHit(f,:)];       
        [dataGroup(subj).hrbytACSbin_FITHit(f,:), dataGroup(subj).hrbytACSbin_RESNORMHit(f,:),dataGroup(subj).hrbytACSbin_RESIDUALHit(f,:),dataGroup(subj).hrbytACSbin_prefPhaseHit(f,:)] = fitCos2RatebyPhase (data2fitHit',X,params);
        dataGroup(subj).hrbytACSbin_FIT_yAllHit(f,:) = dataGroup(subj).hrbytACSbin_FITHit(f,2) + dataGroup(subj).hrbytACSbin_FITHit(f,3).*(cos(X + dataGroup(subj).hrbytACSbin_FITHit(f,1))); % if you want to plot the predicted function
        
        %fit by tACS phase RT
        data2fitRT = [dataGroup(subj).HR_PhaseBinRT(f,:) dataGroup(subj).HR_PhaseBinRT(f,:)];
        [dataGroup(subj).hrbytACSbin_FITRT(f,:), dataGroup(subj).hrbytACSbin_RESNORMRT(f,:),dataGroup(subj).hrbytACSbin_RESIDUALRT(f,:),dataGroup(subj).hrbytACSbin_prefPhaseRT(f,:)] = fitCos2RatebyPhase (data2fitRT',X,params);
        dataGroup(subj).hrbytACSbin_FIT_yAllRT(f,:) = dataGroup(subj).hrbytACSbin_FITRT(f,2) + dataGroup(subj).hrbytACSbin_FITRT(f,3).*(cos(X + dataGroup(subj).hrbytACSbin_FITRT(f,1))); % if you want to plot the predicted function
           
            %create surrogate distribution by shuffling the hit vs miss label
            for ite =1:1000
                surrogGaps = currGaps(randperm(length(currGaps)),5);
               close all
                for bin = 1:nBins
                    if any(binIdx==bin)
                        dataGroup(subj).HR_PhaseBinHitS(f,bin,ite) = mean(surrogGaps(binIdx==bin));
                       % figure,circ_plot(surrogGaps(binIdx==bin))
                    end
                end
                %fit by tACS phase Hit
                data2fitHitS = [dataGroup(subj).HR_PhaseBinHitS(f,:,ite) dataGroup(subj).HR_PhaseBinHitS(f,:,ite)];
                [dataGroup(subj).hrbytACSbin_FITHitS(f,:,ite), dataGroup(subj).hrbytACSbin_RESNORMHitS(f,:,ite),dataGroup(subj).hrbytACSbin_RESIDUALHitS(f,:,ite),dataGroup(subj).hrbytACSbin_prefPhaseHitS(f,:,ite)] = fitCos2RatebyPhase (data2fitHitS',X,params);
            end
            
            %compute z-score relative to surrogate distribution and compute p-value
            
            realvsSurrog_zcore(subj,f) = (dataGroup(subj).hrbytACSbin_FITHit(f,3)-nanmean(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:)))/std(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:));
            realvsSurrog_Pval(subj,f) = mean(abs(dataGroup(subj).hrbytACSbin_FITHit(f,3))<abs(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:)));
            realvsSurrog_percent(subj,f) = (dataGroup(subj).hrbytACSbin_FITHit(f,3)-mean(abs(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:))))./mean(abs(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:)));
      %p = mean(abs(true_value) < abs(shuffled_iterations));
            %compute confidence intervals

SEM = std(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:))/sqrt(length(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:)));               % Standard Error
ts = tinv([0.025  0.975],length(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:))-1);      % T-Score
CIgroup(subj,f,:) = mean(dataGroup(subj).hrbytACSbin_FITHitS(f,3,:)) + ts*SEM;          
        end
        
    
        
    
        %fit hypothetical signals  
        for c =1:size(dataGroup(subj).HR_ImgCosHit,3)
            data2fitTimeHit = [dataGroup(subj).HR_ImgCosHit(f,:,c) dataGroup(subj).HR_ImgCosHit(f,:,c)];  
            [dataGroup(subj).hrbytACSbin_FITTimeHit(f,:,c), dataGroup(subj).hrbytACSbin_RESNORMTimeHit(f,:,c),dataGroup(subj).hrbytACSbin_RESIDUALTimeHit(f,:,c),dataGroup(subj).hrbytACSbin_prefPhaseTimeHit(f,:,c)] = fitCos2RatebyPhase (data2fitTimeHit',XTime,params);
            dataGroup(subj).hrbytACSbin_FIT_yAllTimeHit(f,:,c) = dataGroup(subj).hrbytACSbin_FITTimeHit(f,2,c) + dataGroup(subj).hrbytACSbin_FITTimeHit(f,3,c).*(cos(XTime + dataGroup(subj).hrbytACSbin_FITTimeHit(f,1,c))); % if you want to plot the predicted function
try
            data2fitTimeRT = [dataGroup(subj).HR_ImgCosRT(f,:,c) dataGroup(subj).HR_ImgCosRT(f,:,c)];  
            [dataGroup(subj).hrbytACSbin_FITTimeRT(f,:,c), dataGroup(subj).hrbytACSbin_RESNORMTimeRT(f,:,c),dataGroup(subj).hrbytACSbin_RESIDUALTimeRT(f,:,c),dataGroup(subj).hrbytACSbin_prefPhaseTimeRT(f,:,c)] = fitCos2RatebyPhase (data2fitTimeRT',XTime,params);
            dataGroup(subj).hrbytACSbin_FIT_yAllTimeRT(f,:,c) = dataGroup(subj).hrbytACSbin_FITTimeRT(f,2,c) + dataGroup(subj).hrbytACSbin_FITTimeRT(f,3,c).*(cos(XTime + dataGroup(subj).hrbytACSbin_FITTimeRT(f,1,c))); % if you want to plot the predicted function
catch
    probleanImgRT(subj,c) = 1;
end
        end

    end
     % across frequencies
     for f = 2:5
            realvsSurrog_zcoreAll(subj,f) = (dataGroup(subj).hrbytACSbin_FITHit(f,3)-mean(dataGroup(subj).hrbytACSbin_FITHitS(2:5,3,:),'all'))/std(dataGroup(subj).hrbytACSbin_FITHitS(2:5,3,:),[],'all');
            realvsSurrog_PvalAll(subj,f) = mean(abs(dataGroup(subj).hrbytACSbin_FITHit(f,3))<abs(dataGroup(subj).hrbytACSbin_FITHitS(2:5,3,:)),'all');
            realvsSurrog_percentAll(subj,f) = (dataGroup(subj).hrbytACSbin_FITHit(f,3)-mean(abs(dataGroup(subj).hrbytACSbin_FITHitS(2:5,3,:)),'all'))./mean(abs(dataGroup(subj).hrbytACSbin_FITHitS(2:5,3,:)),'all');
     end
   
    dataHilbert.trCurrGapTrack(:,11) = dataHilbert.tACS1angleGap;
    data2Fit = dataHilbert.trCurrGapTrack;
    data2Fit = data2Fit(~isnan(data2Fit(:,8)),:); %remove nans
    dataGroup(subj).data2Fit = data2Fit;
end

anova_rm(GroupMeanHR)
anova_rm(GroupMeanRT)

%organize some variables
for subj =1:length(dataGroup)
    GroupAmpPhas(subj,:) = dataGroup(subj).hrbytACSbin_FITHit(:,3);
    GroupHRbyPhas_OptPh(subj,:) = dataGroup(subj).hrbytACSbin_prefPhaseHit;
    groupAmp_stimlocked(subj,:,:) = dataGroup(subj).hrbytACSbin_FITTimeHit(:,3,:);
end

%FIGURES FOR PUBLICATION
%Figure 1: Single subject data
iniColorAll = [230/255 0/255 148/255; 186/255 31/255 181/255; 102/255 0/255 161/255; 41/255 5/255 161/255];
upX = linspace(0,2*pi,1000);
figure
for s =1:size(dataGroup,2)
    subplot(4,6,s)
    for currFM = 1:4
        iniColor=iniColorAll(currFM,:);
        plot(X,dataGroup(s).HR_PhaseBinHit(currFM+1,:),'.','Color',iniColor)
        hold on
        yAll(currFM,:) = dataGroup(s).hrbytACSbin_FITHit(currFM+1,2) + dataGroup(s).hrbytACSbin_FITHit(currFM+1,3).*(cos(upX + dataGroup(s).hrbytACSbin_FITHit(currFM+1,1))); % if you want to plot the predicted function
        plot(upX,yAll(currFM,:),'Color',iniColor,'LineStyle','-.')
    end
    %     phases = 0:2*pi/phaseBins:2*pi-2*pi/phaseBins;
    %     meanphases = round(phases+2*pi/phaseBins/2,2);
    xlabel ('Phase')
    ylabel ('Hit rate')
    
    set(gca, 'XTick',X)
    xticklabels({'0' ,'', '', '', '', '2*pi'})
    title (['HR by tACS phase P' num2str(s)])
    %ylim([0.2 1])
end
legend({'0.8 Hz';''; '2';''; '3.2';''; '4.4'})
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',11,'FontWeight','Bold', 'LineWidth', 2);

%Figure 2
%Effect of tACS frequency on HR and reaction time
figure, subplot(1,2,1)
distributionPlotYCC(GroupMeanHR, 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('mean HIT RATE per tACS condition')
xlabel('tACS condition (Hz)')
xticklabels({'Sham' '0.8' '2' '3.2' '4.4'})
ylabel('HIT RATE')
subplot(1,2,2)
distributionPlotYCC(GroupMeanRT, 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('mean REACTION TIME per tACS condition')
xlabel('tACS condition (Hz)')
xticklabels({'Sham' '0.8' '2' '3.2' '4.4'})
ylabel('Reaction Time (s)')

%Figure 3
%amplitude fit plus modulation floor
%Plot mean and standard error of the phase reseted oscillator per condition
%per frequency

KeepbyHR = sum(GroupMeanHR>1,2)==0; %exclude participants with accuracy higher than 0.85 or lower than 0.25 as this would prevent from estimating the amplitude parameter

figure
colorT = [0 1 0; 1 0 0; 0 0 1; 0.3 0.3 0.3; 0 0.6 0.8 ];
xValues = [0.8 2 3.2 4.4]
plot(freqImg,squeeze(mean(groupAmp_stimlocked(KeepbyHR,1,:)))','color',colorT(1,:))
hold on
shadedErrorBar(freqImg,mean(squeeze(groupAmp_stimlocked(KeepbyHR,1,:)))',std(squeeze(groupAmp_stimlocked(KeepbyHR,1,:)))./sqrt(sum(KeepbyHR))');
distributionPlot(GroupAmpPhas(KeepbyHR,2:5),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xValues',freqAll(2:5))%,'xyOri','flipped')
title('mod Amp by tACS Freq')
xlabel('frequency')
set(gca, 'XTick',[0.8 2 3.2 4.4])
xticklabels({'0.8' ,'2', '3.2', '4.4'})
ylabel('amplitude')

%get actual frequency
nbF1 = freqImg==freqAll(2);
nbF2 = freqImg==freqAll(3);
nbF3 = freqImg==freqAll(4);
nbF4 = freqImg==freqAll(5);

% ampTAB = table([mean(groupAmp_stimlocked(KeepbyHR,1,nbF1),3); mean(groupAmp_stimlocked(KeepbyHR,1,nbF2),3);mean(groupAmp_stimlocked(KeepbyHR,1,nbF3),3);mean(groupAmp_stimlocked(KeepbyHR,1,nbF4),3);...
%     GroupAmpPhas(KeepbyHR,2);GroupAmpPhas(KeepbyHR,3);GroupAmpPhas(KeepbyHR,4);GroupAmpPhas(KeepbyHR,5)],[ones(sum(KeepbyHR)*4,1);ones(sum(KeepbyHR)*4,1)+1],[ones(sum(KeepbyHR),1)-0.2;ones(sum(KeepbyHR),1)+1;ones(sum(KeepbyHR),1)+2.2;ones(sum(KeepbyHR),1)+3.4;ones(sum(KeepbyHR),1)-0.2;ones(sum(KeepbyHR),1)+1;ones(sum(KeepbyHR),1)+2.2;ones(sum(KeepbyHR),1)+3.4],repmat(1:sum(KeepbyHR),1,8)',...
%     'VariableNames',{'amptACS','ShamvsReal','FREQ','SUBJ'});
% %ampTAB.FREQ = categorical(ampTAB.FREQ);
% ampTAB.ShamvsReal = categorical(ampTAB.ShamvsReal);
% ampTAB.SUBJ = categorical(ampTAB.SUBJ);
% model1 = fitglme(ampTAB,'amptACS ~ 1 + ShamvsReal + FREQ + (1|SUBJ)','Distribution','normal','FitMethod','Laplace', 'DummyVarCoding', 'effects');
% model2 = fitglme(ampTAB,'amptACS ~ 1 + ShamvsReal*FREQ + (1|SUBJ)','Distribution','normal','FitMethod','Laplace', 'DummyVarCoding', 'effects');
% model3 = fitglme(ampTAB,'amptACS ~ 1 + ShamvsReal*FREQ + (FREQ|SUBJ)','Distribution','normal','FitMethod','Laplace', 'DummyVarCoding', 'effects');
% model4 = fitglme(ampTAB,'amptACS ~ 1 + ShamvsReal*FREQ + (ShamvsReal*FREQ|SUBJ)','Distribution','normal','FitMethod','Laplace', 'DummyVarCoding', 'effects');
% 
% modelsAIC = [model1.ModelCriterion.AIC; model2.ModelCriterion.AIC; model3.ModelCriterion.AIC; model4.ModelCriterion.AIC];
% modelsAIC(:,2)=modelsAIC(:,1)-min(modelsAIC(:,1));
% 
% LR = 2*(model2.LogLikelihood - model1.LogLikelihood); % has a X2 distribution with a df equals to number of constrained parameters, here: 1
% pval = 1 - chi2cdf(LR, 1);

tabletACS = table(mean(groupAmp_stimlocked(KeepbyHR,1,nbF1),3), mean(groupAmp_stimlocked(KeepbyHR,1,nbF2),3),mean(groupAmp_stimlocked(KeepbyHR,1,nbF3),3),mean(groupAmp_stimlocked(KeepbyHR,1,nbF4),3),...
     GroupAmpPhas(KeepbyHR,2),GroupAmpPhas(KeepbyHR,3),GroupAmpPhas(KeepbyHR,4),GroupAmpPhas(KeepbyHR,5),...
'VariableNames',{'sham08','sham2','sham32','sham44','tACS08','tACS2','tACS32','tACS44'});
withinDesign = table(([1 2 3 4 1 2 3 4])',([1 1 1 1 2 2 2 2])','VariableNames',{'Freq','Cond'});
withinDesign.Freq = categorical(withinDesign.Freq);
withinDesign.Cond = categorical(withinDesign.Cond);
%rmfitting - table - testing groups - grouping
myRM = fitrm(tabletACS,'sham08-tACS44 ~ 1','withinDesign',withinDesign);
 
[b1_rm,A_rm,C_rm,D_rm] = ranova(myRM,'WithinModel','Freq*Cond');
Multicompare_tableTACS = multcompare(myRM,'Freq');
Multicompare_tableMD = multcompare(myRM,'Cond');

SEM_sham08 = std(mean(groupAmp_stimlocked(KeepbyHR,1,nbF1),3))/sqrt(length(mean(groupAmp_stimlocked(KeepbyHR,1,nbF1),3)));               % Standard Error
ts_sham08 = tinv([0.025  0.975],length(mean(groupAmp_stimlocked(KeepbyHR,1,nbF1),3))-1);      % T-Score
CI_sham08 = mean(mean(groupAmp_stimlocked(KeepbyHR,1,nbF1),3)) + ts_sham08*SEM_sham08;                      % Confidence Intervals

SEM_sham2 = std(mean(groupAmp_stimlocked(KeepbyHR,1,nbF2),3))/sqrt(length(mean(groupAmp_stimlocked(KeepbyHR,1,nbF2),3)));               % Standard Error
ts_sham2 = tinv([0.025  0.975],length(mean(groupAmp_stimlocked(KeepbyHR,1,nbF2),3))-1);      % T-Score
CI_sham2 = mean(mean(groupAmp_stimlocked(KeepbyHR,1,nbF2),3)) + ts_sham2*SEM_sham2;                      % Confidence Intervals

SEM_sham32 = std(mean(groupAmp_stimlocked(KeepbyHR,1,nbF3),3))/sqrt(length(mean(groupAmp_stimlocked(KeepbyHR,1,nbF3),3)));               % Standard Error
ts_sham32 = tinv([0.025  0.975],length(mean(groupAmp_stimlocked(KeepbyHR,1,nbF3),3))-1);      % T-Score
CI_sham32 = mean(mean(groupAmp_stimlocked(KeepbyHR,1,nbF3),3)) + ts_sham32*SEM_sham32;                      % Confidence Intervals

SEM_sham44 = std(mean(groupAmp_stimlocked(KeepbyHR,1,nbF4),3))/sqrt(length(mean(groupAmp_stimlocked(KeepbyHR,1,nbF4),3)));               % Standard Error
ts_sham44 = tinv([0.025  0.975],length(mean(groupAmp_stimlocked(KeepbyHR,1,nbF4),3))-1);      % T-Score
CI_sham44 = mean(mean(groupAmp_stimlocked(KeepbyHR,1,nbF4),3)) + ts_sham44*SEM_sham44;                      % Confidence Intervals

sigPart1 = find(GroupAmpPhas(KeepbyHR,2)>CI_sham08(2));
sigPart2 = find(GroupAmpPhas(KeepbyHR,3)>CI_sham2(2));
sigPart3 = find(GroupAmpPhas(KeepbyHR,4)>CI_sham32(2));
sigPart4 = find(GroupAmpPhas(KeepbyHR,5)>CI_sham44(2));

sigPart1neg = find(GroupAmpPhas(KeepbyHR,2)<CI_sham08(1));
sigPart2neg = find(GroupAmpPhas(KeepbyHR,3)<CI_sham2(1));
sigPart3neg = find(GroupAmpPhas(KeepbyHR,4)<CI_sham32(1));
sigPart4neg = find(GroupAmpPhas(KeepbyHR,5)<CI_sham44(1));

%plot significance per frequency as a matrix for visualization
sigPosPart = zeros(sum(KeepbyHR),4);
signegPart = zeros(sum(KeepbyHR),4);

for s =1:sum(KeepbyHR)
    if sum(sigPart1==s)
        sigPosPart(s,1) = 1;
    end
    if sum(sigPart2==s)
        sigPosPart(s,2) = 1;
    end
    if sum(sigPart3==s)
        sigPosPart(s,3) = 1;
    end
    if sum(sigPart4==s)
        sigPosPart(s,4) = 1;
    end
    
    %negative
     if sum(sigPart1neg==s)
        signegPart(s,1) = 1;
    end
    if sum(sigPart2neg==s)
        signegPart(s,2) = 1;
    end
    if sum(sigPart3neg==s)
        signegPart(s,3) = 1;
    end
    if sum(sigPart4neg==s)
        signegPart(s,4) = 1;
    end
end
figure
for s =1:length(sigPosPart)
subplot(4,6,s)
imagesc([sigPosPart(s,:);signegPart(s,:)])

end

figure
subplot(1,2,1)
imagesc(sigPosPart)

subplot(1,2,2)
imagesc(signegPart)

%figure 4
figure
xValues2 = xValues +0.4;
distributionPlot([mean(groupAmp_stimlocked(KeepbyHR,1,nbF1),3) GroupAmpPhas(KeepbyHR,2) mean(groupAmp_stimlocked(KeepbyHR,1,nbF2),3) GroupAmpPhas(KeepbyHR,3) mean(groupAmp_stimlocked(KeepbyHR,1,nbF3),3) GroupAmpPhas(KeepbyHR,4) mean(groupAmp_stimlocked(KeepbyHR,1,nbF4),3) GroupAmpPhas(KeepbyHR,5)],'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('mod Amp by Freq')
xlabel('frequency')
set(gca, 'XTick',[0.8 2 3.2 4.4])
xticklabels({'0.8' ,'2', '3.2', '4.4'})
ylabel('amplitude')

ttest(mean(groupAmp_stimlocked(KeepbyHR,1,nbF1),3), GroupAmpPhas(KeepbyHR,2))
ttest(mean(groupAmp_stimlocked(KeepbyHR,1,nbF2),3), GroupAmpPhas(KeepbyHR,3))
ttest(mean(groupAmp_stimlocked(KeepbyHR,1,nbF3),3), GroupAmpPhas(KeepbyHR,4))
ttest(mean(groupAmp_stimlocked(KeepbyHR,1,nbF4),3), GroupAmpPhas(KeepbyHR,5))


%figure 5, optimal tACS phase
figure
subplot(1,4,1)
circ_plot(GroupHRbyPhas_OptPh(KeepbyHR,2),'pretty')
title('0.8')
subplot(1,4,2)
circ_plot(GroupHRbyPhas_OptPh(KeepbyHR,3),'pretty')
title('2')
subplot(1,4,3)
circ_plot(GroupHRbyPhas_OptPh(KeepbyHR,4),'pretty')
title('3.2')
subplot(1,4,4)
circ_plot(GroupHRbyPhas_OptPh(KeepbyHR,5),'pretty')
title('Pref tACS phase by Freq 4.4')
legend

[p1,z1] = circ_rtest (GroupHRbyPhas_OptPh(:,2));
[p2,z2] = circ_rtest (GroupHRbyPhas_OptPh(:,3));
[p3,z3] = circ_rtest (GroupHRbyPhas_OptPh(:,4));
[p4,z4] = circ_rtest (GroupHRbyPhas_OptPh(:,5));

%compute percent change

percentChange = [(tabletACS.tACS08-tabletACS.sham08)./tabletACS.sham08,...
    (tabletACS.tACS2-tabletACS.sham2)./tabletACS.sham2,...
    (tabletACS.tACS32-tabletACS.sham32)./tabletACS.sham32,...
    (tabletACS.tACS44-tabletACS.sham44)./tabletACS.sham44];
    
figure
xValues2 = xValues +0.4;
distributionPlot(percentChange,'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('mod Amp by Freq')
xlabel('frequency')
set(gca, 'XTick',[0.8 2 3.2 4.4])
xticklabels({'0.8' ,'2', '3.2', '4.4'})
ylabel('amplitude')

ttest(percentChange(:,1))
ttest(percentChange(:,2))
ttest(percentChange(:,3))
ttest(percentChange(:,4))

figure, subplot(1,2,1)
distributionPlotYCC(realvsSurrog_zcoreAll, 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('z-score')
xlabel('tACS condition (Hz)')
xticklabels({'Sham' '0.8' '2' '3.2' '4.4'})
subplot(1,2,2)
distributionPlotYCC(realvsSurrog_zcore, 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('zscore')
xlabel('tACS condition (Hz)')
xticklabels({'Sham' '0.8' '2' '3.2' '4.4'})

[h,p,ci,stats1] = ttest(realvsSurrog_zcoreAll(:,2))
[h,p,ci,stats1] = ttest(realvsSurrog_zcoreAll(:,3))
[h,p,ci,stats1] = ttest(realvsSurrog_zcoreAll(:,4))
[h,p,ci,stats1] = ttest(realvsSurrog_zcoreAll(:,5))

[h1,p1,ci1,stats11] = ttest(realvsSurrog_zcore(:,2))
[h2,p2,ci2,stats12] = ttest(realvsSurrog_zcore(:,3))
[h3,p3,ci3,stats13] = ttest(realvsSurrog_zcore(:,4))
[h4,p4,ci4,stats14] = ttest(realvsSurrog_zcore(:,5))


save (fullfile(OUTdir,'groupData_tACSonlyNoise_InManuscript_Jan2025'), '-v7.3')

