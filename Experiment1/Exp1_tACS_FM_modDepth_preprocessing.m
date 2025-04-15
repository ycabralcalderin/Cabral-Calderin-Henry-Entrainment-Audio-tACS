% Pipeline for analysing FM-stimulus+tACS data
% this script reads the eeg raw data file and prepares a matlab variable
% including for each subject each gap response together with FM and tACS
% variables
% Cabral-Calderin May 2022. Last modified January. 2025

% Experiment 1- Cabral-Calderin & Henry, Plos Biology 2025
% Additional Toolboxes needed: CircStat2012a, fieldtrip-20201128

clear
close all
clc

%DEFINE DIRECTORIES
MAINdir  = uigetdir (); %main path to raw data, e.g., '/Users/yuranny.cabral/Documents/gapDetectionProject/tACS/DATA'
BEHdir   = fullfile(MAINdir, 'rawdata/Behavioral'); %raw beh data
EEGdir   = fullfile(MAINdir, 'rawdata/tACS'); %raw eeg data (for tACS signal)
OUTdir   = fullfile(MAINdir, 'ANA'); %main output directory
addpath('/CircStat2012a')
addpath('/fieldtrip-20201128/');
addpath(genpath('/Scripts/tACS')) %own scripts

%define subject and protocol info % Subject ID, FM depth S1, FMdepth S2, FM
%depth S3, number of sessions
SUBJlist = {'1' 11 39 [] 2; '2' 11 39 [] 2; '3' 39 11 [] 2; '4' 39 11 [] 2; '5' 39 11 [] 2;
            '6' 39 11 [] 2; '7' 11 39 [] 2; '8' 39 [] [] 1; '9' 39 11 [] 2;'10' 11 39 [] 2;
           '11' 11 39 [] 2;'12' 39 [] [] 1;'13' 11 39 [] 2;'14' 11 39 [] 2;'15' 11,39 [] 2;
           '16' 39 11 [] 2;'17' 11 39 [] 2;'18' 39 11 [] 2;'19' 11 39 [] 2;'20' 39 11 [] 2;
           '21' 11 [] [] 1;'22' 11 39 [] 2;'23' 11 39 [] 2;'24' 39 11 [] 2;'25' 11 39 [] 2;
           '26' 39 [] [] 1;'27' 39 11 [] 2;'28' 39 11 [] 2;'29' 11 39 [] 2;'30' 11 39 [] 2;
           '31' 11 39 [] 2;'32' 39 11 [] 2;'33' 39 11 [] 2;'34' 11 39 [] 2;};

ft_defaults %initialize fieldtrip defaults

%Initialize some variables
dataGroupallGaps = struct([]);
checkSubject     = zeros(length(SUBJlist),3);
gapSizeGroup     = struc([]);

for subj = 1:length(SUBJlist)%loop across subjects
    close all
    currSUBJ = SUBJlist{subj,1}; % current subject
    fprintf('%s%s%/nn','... working on subject ', currSUBJ)
    
    %get all files for the current participant
    files = dir([fullfile(EEGdir,currSUBJ) '_*.eeg']);
    
    %loop through sessions
    for Session = 1:SUBJlist{subj,5}%length(files)
        
        sessionID = str2double(files(Session).name(end-4));
        if strcmp(currSUBJ,'5')&& strcmp(files(Session).name,'5_tACS1_block1-4.eeg')
            sessionID = 1;
        elseif strcmp(currSUBJ,'29')&& Session == 2
            sessionID = str2double(files(Session).name(end-6));
        end
        
        %load stim file to check that the modulation depth is correct
        if strcmp(currSUBJ,'11')&& Session == 1
            load(fullfile(BEHdir,currSUBJ,num2str(Session), '1/MAIN/STIM',[currSUBJ '_Stim_Main_orig.mat']));
        elseif strcmp(currSUBJ,'25')&& Session == 1 %there were two folders, one for the first 3 blocks and one for the last 2
            load(fullfile(BEHdir,currSUBJ,num2str(Session),'block 4 and 5/STIM/MKY03_Stim_Main_orig.mat'));
        else
            load(fullfile(BEHdir,currSUBJ,num2str(Session), 'MAIN/STIM',[currSUBJ '_Stim_Main_orig.mat']));
        end
        
        if isfield(Stim(1).stiminfo,'ModDepth')
            if SUBJlist{subj,Session+1}~=Stim(1).stiminfo.ModDepth
                error 'WARNING: Modulation depth in SUBJlist doesnt match with data'
            end
        else
        end
        if Session ~= sessionID
            error 'CHECK SESSION AND SESSION ID DO NOT MATCH'
        end
        curroDIR = fullfile(OUTdir,currSUBJ); %subject dir for saving data
        
        if exist(curroDIR,'dir') %create dir if it doesn't exist
        else
            mkdir(curroDIR)
        end
        
        %get gap size
        gapSizeGroup(subj,Session) = Stim(1).stiminfo.gapDur;
        
        %WORK on the .eeg files to get the tACS signal and triggers and define
        %tACS at stimulus time
        
        if strcmp(currSUBJ,'5')&& strcmp(files(Session).name,'5_tACS1_block1-4.eeg')
            sessionEEG = '1_block1-4'; %for this participant, only blocks 1-4 were collected in the EEG. The last block was Sham. Load matlab files to get this data
            
            load(fullfile(BEHdir,currSUBJ,num2str(Session), 'MAIN',[currSUBJ '_Block_5_ongoingTask.mat']));
            
            %get STIM Phase and detect for each trial
            ShamData = [];
            for tr = 129:160
                for g = 1:length(response(tr).target)
                    if isempty(response(tr).target(g).r)
                        ShamData = [ShamData; response(tr).target(g).bin 0];
                    else
                        ShamData = [ShamData; response(tr).target(g).bin 1];
                    end
                end
            end
        else
            sessionEEG = sessionID;
        end
        if strcmp(currSUBJ,'29')&& Session ==2 %this subject had two files in the second session
            
            [dataHilbert, lphpData_Stim,dataCutGaps]    = get_tACSphaseStimpreV2(fullfile(EEGdir,currSUBJ),'2_1');
            [dataHilbert2, lphpData_Stim2,dataCutGaps2] = get_tACSphaseStimpreV2(fullfile(EEGdir,currSUBJ),'2_2');
            
            %concatenate data from both files
            dataHilbert2.trCurrStimTrack(:,1:2) = dataHilbert2.trCurrStimTrack(:,1:2)+length(dataHilbert.tACS1); %add length of first file to have consecutive numbers
            dataHilbert2.trCurrGapTrack(:,1:2)  = dataHilbert2.trCurrGapTrack(:,1:2)+length(dataHilbert.tACS1);
            F = fieldnames(dataHilbert); % Return an alphanumerically sorted listing
            for ix=1:numel(F)
                if size(dataHilbert.(F{ix}),1)>size(dataHilbert.(F{ix}),2)
                    dataHilbert.(F{ix}) = [dataHilbert.(F{ix}); dataHilbert2.(F{ix})];
                else
                    dataHilbert.(F{ix}) = [dataHilbert.(F{ix}) dataHilbert2.(F{ix})];
                end
            end
            %add last sample info from file 1 to avoid sample issues when
            %concatenating
            lphpData_Stim2.sampleinfo = lphpData_Stim2.sampleinfo +length(dataHilbert.tACS1);
            lphpData_Stim             = ft_appenddata([], lphpData_Stim,lphpData_Stim2);
            
            dataCutGaps2.sampleinfo = dataCutGaps2.sampleinfo + length(dataHilbert.tACS1);
            dataCutGaps             = ft_appenddata([], dataCutGaps,dataCutGaps2);
        else
            [dataHilbert, lphpData_Stim,dataCutGaps] = get_tACSphaseStimpreV2(fullfile(EEGdir,currSUBJ),sessionEEG);
        end
        if length(dataHilbert.trCurrStimTrack)>160 %just double check in case there are more blocks than expected
            fprintf ('WARNING, more than 5 blocks detected')
            disp(checkSubject(subj,Session))
        end
        
        %concatenate trials
        alltrials = cat(1,dataCutGaps.trial{:});
        channel1  = 1:2:size(alltrials,1);
        tdiff     = (max(alltrials(channel1,:),[],2)- min(alltrials(channel1,:),[],2))'; %get range per trial: max amp-min
        ThrAmp    = min(tdiff)+((median(tdiff)-min(tdiff))./3); %threshold to decide if it was a stim or a sham trial
        
        BlockLim       = 1:136:length(dataCutGaps.trial);
        if strcmp(currSUBJ,'5')&& strcmp(files(Session).name,'5_tACS1_block1-4.eeg')
            blockmeanrange = [mean(tdiff(BlockLim(1):BlockLim(2)-1)) mean(tdiff(BlockLim(2):BlockLim(3)-1)) mean(tdiff(BlockLim(3):BlockLim(4)-1)) mean(tdiff(BlockLim(4):length(dataCutGaps.trial)))];
            firstSham      = length(dataHilbert.trCurrGapTrack(:,5))+1;
            lastSham       = firstSham +135; %each block has about 136 gaps
            tdiff          = [tdiff zeros(1,136)];%add zeros for the sham gap trials, which were not recorded in the EEG
            NewStimTr      = tdiff>=ThrAmp;
            NewShamTr1     = tdiff<ThrAmp;
            NewShamTr      = zeros(1,length(NewShamTr1));
            NewShamTr(firstSham:lastSham) = NewShamTr1(firstSham:lastSham); %take sham trials only within the sham block
            NewStimTr(firstSham:lastSham) = 0; %replace possible sham trials wrongly considered as stim
        else
            if strcmp(currSUBJ,'11')&& Session ==1 ||strcmp(currSUBJ,'31')&& Session ==1  %this participant had one extra block
                blockmeanrange = [mean(tdiff(BlockLim(1):BlockLim(2)-1)) mean(tdiff(BlockLim(2):BlockLim(3)-1)) mean(tdiff(BlockLim(3):BlockLim(4)-1)) mean(tdiff(BlockLim(4):BlockLim(5)-1)) mean(tdiff(BlockLim(5):BlockLim(6)-1)) mean(tdiff(BlockLim(6):length(dataCutGaps.trial)))];
            elseif strcmp(currSUBJ,'5')&& Session ==1  %this participant was missing one block
                blockmeanrange = [mean(tdiff(BlockLim(1):BlockLim(2)-1)) mean(tdiff(BlockLim(2):BlockLim(3)-1)) mean(tdiff(BlockLim(3):BlockLim(4)-1)) mean(tdiff(BlockLim(4):length(dataCutGaps.trial)))];
            else
                blockmeanrange = [mean(tdiff(BlockLim(1):BlockLim(2)-1)) mean(tdiff(BlockLim(2):BlockLim(3)-1)) mean(tdiff(BlockLim(3):BlockLim(4)-1)) mean(tdiff(BlockLim(4):BlockLim(5)-1)) mean(tdiff(BlockLim(5):length(dataCutGaps.trial)))];
            end
            [minrange,idx] = min(blockmeanrange);%get block with smallest range as the sham
            
            %firstSham = str2double(inputdlg('first sham trial')); %select the first sham trial based on the channel data
            if strcmp(currSUBJ,'11')&& Session ==1 %there were two sham trials
                idx = find(blockmeanrange<50);%get block with smallest range as the sham
            end
            firstSham  = BlockLim(idx);
            lastSham   = firstSham +135; %each block has about 136 gaps
            NewStimTr  = tdiff>=ThrAmp;
            NewShamTr1 = tdiff<ThrAmp;
            NewShamTr  = zeros(1,length(NewShamTr1));
            if strcmp(currSUBJ,'30')&& Session ==2 %a few trials were missing at the end
                lastSham = 672;
            elseif strcmp(currSUBJ,'1')&& Session ==2 %a few trials were missing at the end
                lastSham = 675;
            end
            for n=1:length(firstSham)
                NewShamTr(firstSham(n):lastSham(n)) = NewShamTr1(firstSham(n):lastSham(n)); %take sham trials only within the sham block
                NewStimTr(firstSham(n):lastSham(n)) = 0; %replace possible sham trials wrongly considered as stim
            end
            
        end
        if strcmp(currSUBJ,'11')&& Session ==1 %there were two sham trials
        else
            if sum(NewShamTr)<122 || sum(NewShamTr)>175
                disp(sum(NewShamTr))
                error 'Sham trials were not correctly identified'
            end
        end
        disp(firstSham)
        
        %do some plotting for visualizing
        
        figure, scatter(1:length(tdiff),tdiff)
        hold on
        plot([1 length(tdiff)],[ThrAmp ThrAmp])
        scatter(find(NewShamTr==1),tdiff(NewShamTr==1)-50,'*g')
        scatter(find(NewStimTr==1),tdiff(NewStimTr==1),'b')
        %pause(2)
        save (fullfile(curroDIR,'preproctACS'),'dataHilbert','lphpData_Stim','NewShamTr','NewStimTr','-v7.3')
        
        AUDbin = 0:2*pi/8:2*pi; %possible starting phases for the audio stim
        
        dataHilbert.trCurrStimTrack(:,6) = AUDbin(dataHilbert.trCurrStimTrack(:,5)); %add phase in radians
        
        tACS_AUDdist = circ_dist(wrapTo2Pi(dataHilbert.tACS1angleStim)',dataHilbert.trCurrStimTrack(:,6));
        
        %reconstruct signal
        f    = 2;
        gapAUDphase = nan(1,length(dataHilbert.trCurrGapTrack));
        AUDTr       = struct([]);
        gapAUDtACSlag = nan(1,length(dataHilbert.trCurrGapTrack));
        
        for Tr=1:length(dataHilbert.trCurrStimTrack)
            
            AUDTr(Tr).Gap = find(dataHilbert.trCurrGapTrack(:,1)>dataHilbert.trCurrStimTrack(Tr,1)&dataHilbert.trCurrGapTrack(:,1)<dataHilbert.trCurrStimTrack(Tr,2));
            gapAUDtACSlag(AUDTr(Tr).Gap) = tACS_AUDdist(Tr); % phase lag between tACS and the audio signal
            
            %get stim-related phase
            t        = dataHilbert.trCurrStimTrack(Tr,1)-dataHilbert.trCurrStimTrack(Tr,3):dataHilbert.trCurrStimTrack(Tr,2);
            t2       = 0:1/5000:(t(end)-t(1))/5000;
            startPha = dataHilbert.trCurrStimTrack(Tr,6);
            AUDPha   = wrapTo2Pi(2*pi*f*t2+startPha);
            yCos     = cos(2*pi*f*t2+startPha);
            
            
            figure,
            plot(yCos)
            hold on
            plot(AUDPha)
            %gap start within the sound
            gapstart = (dataHilbert.trCurrGapTrack(AUDTr(Tr).Gap,1)-dataHilbert.trCurrGapTrack(AUDTr(Tr).Gap,3))-(dataHilbert.trCurrStimTrack(Tr,1)-dataHilbert.trCurrStimTrack(Tr,3));
            gapAUDphase(AUDTr(Tr).Gap) = AUDPha(gapstart);
            
            hold off
        end
        if strcmp(currSUBJ,'5')&& strcmp(files(Session).name,'5_tACS1_block1-4.eeg')
            dataHilbert.trCurrGapTrack(length(dataHilbert.trCurrGapTrack(:,4:5))+1:length(dataHilbert.trCurrGapTrack(:,4:5))+length(ShamData),4:5) = ShamData;
            dataHilbert.tACS1angleGap = [dataHilbert.tACS1angleGap nan(1,length(ShamData))];
            X         = 0:2*pi/9:2*pi-2*pi/9; %this is my phase vector
            X         = X+2*pi/9/2-0.1214;% the 0.1214 is to correct for the difference between finding the phase at the middle of the gap vs at the beginning as is done using the EEG data
            
            gapAUDphase = [gapAUDphase X(ShamData(:,1))];
            gapAUDtACSlag = [gapAUDtACSlag nan(1,length(ShamData))];
        end
        %put together data for logistic regression later on
        dataGroupallGaps(subj,Session).firstShamTr    = firstSham;
        dataGroupallGaps(subj,Session).gapDetect      = dataHilbert.trCurrGapTrack(:,5);
        dataGroupallGaps(subj,Session).gaptACSangle   = dataHilbert.tACS1angleGap;
        dataGroupallGaps(subj,Session).gapAUDangle    = gapAUDphase;
        dataGroupallGaps(subj,Session).gapAUD_tACSlag = gapAUDtACSlag;
        dataGroupallGaps(subj,Session).gapStimcond    = zeros(length(gapAUDphase),1); %zero mean undefined
        dataGroupallGaps(subj,Session).gapStimcond(NewStimTr) = 1; %1 means real stim
        dataGroupallGaps(subj,Session).gapStimcond(logical(NewShamTr)) = 2; %2 means sham
        dataGroupallGaps(subj,Session).ModDepth       = SUBJlist{subj,Session+1};
        dataGroupallGaps(subj,Session).SessionNr      = Session;
        dataGroupallGaps(subj,Session).FA             = dataHilbert.FA;
        try
            dataGroupallGaps(subj,Session).FAtACSlag      = tACS_AUDdist(dataHilbert.FAstim);
        catch
            dataGroupallGaps(subj,Session).FAtACSlag = 0;
        end
        dataGroupallGaps(subj,Session).resp2Late      = dataHilbert.resp2Late;
        dataGroupallGaps(subj,Session).resp2LateRT    = dataHilbert.resp2LateRT;
        try
            dataGroupallGaps(subj,Session).resp2LatetACSlag = tACS_AUDdist(dataHilbert.resp2Latestim);
        catch
            dataGroupallGaps(subj,Session).resp2LatetACSlag = 0;
        end
        %Change session order to have it according to the Mod Depth
        if dataGroupallGaps(subj,Session).ModDepth == 11
            allDataModDepth(subj,1) = dataGroupallGaps(subj,Session);
        elseif dataGroupallGaps(subj,Session).ModDepth == 39
            
            allDataModDepth(subj,2) = dataGroupallGaps(subj,Session);
            
        end
        
        close all
        
    end
end

save (fullfile(OUTdir,'dataGroupallGaps_tACSModDepth_Jan2025'),'dataGroupallGaps','allDataModDepth','gapSizeGroup','-v7.3')


