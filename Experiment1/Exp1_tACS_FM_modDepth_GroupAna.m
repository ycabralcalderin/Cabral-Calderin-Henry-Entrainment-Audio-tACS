%pipeline for analysing FM-stimulus+tACS data
clear
close all
clc

%DEFINE DIRECTORIES
OUTdir   = fullfile('/mnt/beegfs/users/yuranny.cabral/FMtACSModDepth/ANA'); %main output directory
cd(OUTdir)

%define some general variables
NrGapBins     = 9; %number of FM stimulus phase bins for gap presentation
NrAudtACSlag  = 6; %tACS-Audio phase lag. 6 chosen based on previous tACS study

%load Group data with the info about detection FM phase and tACS phase.
%group variables were created in "ANApipeline_tACSFMStim_PrepareGroupVariable_ModDepth.m" 
load(fullfile(OUTdir,'dataGroupallGaps_tACSModDepth_March2024'),'dataGroupallGaps','allDataModDepth','gapSizeGroup')
%define subject and protocol info
SUBJlist = {'AAH07' 11 39 [] 2;'ABI24' 11 39 [] 2;'ADA23' 39 11 [] 2;'ALC04' 39 11 [] 2;'ASE25' 39 11 [] 2;
    'ATE26' 39 11 [] 2;'BAO30' 11 39 [] 2;'BSD07' 39 [] [] 1;'BWC24' 39 11 [] 2;
    'CSA27' 11 39 [] 2;'EGA23' 11 39 [] 2;'EHA04' 39 [] [] 1;'FAL19' 11 39 [] 2;'IAN18' 11 39 [] 2;
    'JBL26' 11,39 [] 2;'JWE28' 39 11 [] 2;'KBA19' 11 39 [] 2;'KWA10' 39 11 [] 2;
    'LCS16' 11 39 [] 2;'LLY08' 39 11 [] 2;'LPH19' 11 [] [] 1;'MKA11' 11 39 [] 2; %subject ID
    'MKA23' 11 39 [] 2;'MKS13' 39 11 [] 2;'MKZ03' 11 39 [] 2;'MPA12' 39 [] [] 1;'MRN19' 39 11 [] 2;
    'NYL03' 39 11 [] 2;'RIA10' 11 39 [] 2;'RRN10' 11 39 [] 2;'SKA14' 11 39 [] 2;
    'SSA10' 39 11 [] 2;'STK21' 39 11 [] 2;'YEN03' 11 39 [] 2;};

%code session order for analysis
depthOrd = zeros(length(SUBJlist),1);
for subj=1:length(SUBJlist)%loop across subjects
    if SUBJlist{subj,2}==11
        depthOrd (subj)=1;
    else
depthOrd (subj)=2;
    end
end

%parameters for binning and for fitting the cosine
phasebinGap = 0:2*pi/NrGapBins:2*pi;
lagBins     = 0:2*pi/NrAudtACSlag:2*pi; %0:2*pi/4:2*pi; %this is the phase in the center of the bin
maxSeparation = 2*pi/NrAudtACSlag/2; %this is  the half width of a bin for not overlapping
lagBins1 =  wrapTo2Pi(lagBins-maxSeparation);
lagBins2 =  wrapTo2Pi(lagBins+maxSeparation);
NrIte = 1000;
%fit Cosine function
X         = 0:2*pi/NrGapBins:2*pi-2*pi/NrGapBins; %this is my phase vector
X         = X+2*pi/NrGapBins/2;
lag       = pi;
intercept = 0.5;
amp       = 0.5;
params    = [lag intercept amp];


%PLOTTING TO VISUALIZE THE PHASE LAG BINS
figure,
for ll = 1:length(lagBins)-1
    subplot(2,3,ll)
    circ_plot([lagBins1(ll) lagBins(ll) lagBins2(ll)]','pretty')
end
                        
for subj = 1:length(SUBJlist)%loop across subjects
    %   close all
    currSUBJ = SUBJlist{subj,1}; % current subject
    fprintf('%s%s\n','... working on subject ', currSUBJ)
    
    NrSess         = SUBJlist{subj,5};
    gapDetect      = [];
    gaptACSangle   = [];
    gapAUDangle    = [];
    gapAUD_tACSlag = [];
    gapStimcond    = [];
    
    for Session = 1:NrSess
        if ~isempty(allDataModDepth(subj,Session).gapDetect)
          
            gapDetect      = allDataModDepth(subj,Session).gapDetect;
            gaptACSangle   = allDataModDepth(subj,Session).gaptACSangle';
            gapAUDangle    = allDataModDepth(subj,Session).gapAUDangle';
            gapAUD_tACSlag = allDataModDepth(subj,Session).gapAUD_tACSlag';
            gapStimcond    = allDataModDepth(subj,Session).gapStimcond;
            
            %first separate sham from real
            %main
            tACS1angleGapMain  = wrapTo2Pi(gaptACSangle(gapStimcond==1)); %tACS phase at gap onset
            trCurrGapTrackMain = gapDetect(gapStimcond==1); %gap info from eeg triggers
            gapAUDtACSlagMain  = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==1)); %phase lag between tACS and audio signal
            gapAUDphaseMain    = gapAUDangle(gapStimcond==1); %FM phase at gap onset
            %sham
            tACS1angleGapsham  = wrapTo2Pi(gaptACSangle(gapStimcond==2)); %tACS phase at gap onset
            trCurrGapTracksham = gapDetect(gapStimcond==2); %gap info from eeg triggers
            gapAUDtACSlagsham  = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==2)); %phase lag between tACS and audio signal
            gapAUDphasesham    = gapAUDangle(gapStimcond==2); %FM phase at gap onset
            
            for gapBin = 1:length(phasebinGap)-1
                currGapbintACS  = find(tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))'; %gaps in current tACS phase bin
                currGapbinStim  = find(gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))'; %gaps im current AUD bin
                                length(currGapbintACS)
                                length(currGapbinStim)
                for b=1:length(lagBins)-1 %separate by tACS-AUD phase lag
                    if lagBins1(b)>lagBins2(b)
                        currGaptACSbinSlag = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&...
                            (gapAUDtACSlagMain>=lagBins1(b) | gapAUDtACSlagMain<lagBins2(b)))';
                        currGapStimbinSlag = find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                            (gapAUDtACSlagMain>=lagBins1(b) | gapAUDtACSlagMain<lagBins2(b)))';
                    else
                        currGaptACSbinSlag = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&...
                            (gapAUDtACSlagMain>=lagBins1(b)&gapAUDtACSlagMain<lagBins2(b)))';
                        currGapStimbinSlag =find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                            (gapAUDtACSlagMain>=lagBins1(b)& gapAUDtACSlagMain<lagBins2(b)))';
                    end
                    dataGroupPool(subj,Session).hrbytACSbinStimLagMainNr(b,gapBin) = length(currGaptACSbinSlag);
                    dataGroupPool(subj,Session).hrbyStimbintACSLagMainNr(b,gapBin) = length(currGapStimbinSlag);
                    dataGroupPool(subj,Session).hrbytACSbinStimLagMain(b,gapBin) = mean(trCurrGapTrackMain(currGaptACSbinSlag));
                    dataGroupPool(subj,Session).hrbyStimbintACSLagMain(b,gapBin) = mean(trCurrGapTrackMain(currGapStimbinSlag));
                    
                                        figure, subplot(2,2,1), circ_plot(tACS1angleGapMain(currGaptACSbinSlag))
                                        title(['TACS gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
                                        subplot(2,2,2), circ_plot(gapAUDtACSlagMain(currGaptACSbinSlag))
                                        title('lag')
                                        subplot(2,2,3), circ_plot(gapAUDphaseMain(currGapStimbinSlag))
                                        title(['AUDIO gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
                                        subplot(2,2,4), circ_plot(gapAUDtACSlagMain(currGapStimbinSlag))
                                        title('lag')
                end
                %separate by both FM and tACS signal
                for gapBin2 = 1:length(phasebinGap)-1
                    currGapbinAT = find((tACS1angleGapMain>=phasebinGap(gapBin)&tACS1angleGapMain<phasebinGap(gapBin+1))&...
                        (gapAUDphaseMain>=phasebinGap(gapBin2)&gapAUDphaseMain<phasebinGap(gapBin2+1)))';
                    
                    %  length(currGapbinAT)
                    dataGroupPool(subj,Session).hrbyStimbytACSbinMain(gapBin2,gapBin) = mean(trCurrGapTrackMain(currGapbinAT));
                    dataGroupPool(subj,Session).hrbyStimbytACSbinMainNr(gapBin2,gapBin) = length(currGapbinAT);
                end
                dataGroupPool(subj,Session).hrbyStimbinMain(gapBin) = mean(trCurrGapTrackMain(currGapbinStim));
                dataGroupPool(subj,Session).hrbytACSbinMain(gapBin) = mean(trCurrGapTrackMain(currGapbintACS));
            end
            
            for gapBin = 1:length(phasebinGap)-1
                
                currGapbintACSsham = find(tACS1angleGapsham>=phasebinGap(gapBin)&tACS1angleGapsham<phasebinGap(gapBin+1))';
                currGapbinStimsham = find(gapAUDphasesham>=phasebinGap(gapBin)&gapAUDphasesham<phasebinGap(gapBin+1))';
                     dataGroupPool(subj,Session).hrbyStimbinshamNr(gapBin) = length(currGapbintACSsham);
                     dataGroupPool(subj,Session).hrbytACSbinshamNr(gapBin) = length(currGapbinStimsham);
                dataGroupPool(subj,Session).hrbyStimbinsham(gapBin) = mean(trCurrGapTracksham(currGapbinStimsham));
                dataGroupPool(subj,Session).hrbytACSbinsham(gapBin) = mean(trCurrGapTracksham(currGapbintACSsham));
                
                %separate by both FM and tACS signal
                for gapBin2 = 1:length(phasebinGap)-1
                    currGapbinATs = find((tACS1angleGapsham>=phasebinGap(gapBin)&tACS1angleGapsham<phasebinGap(gapBin+1))&...
                        (gapAUDphasesham>=phasebinGap(gapBin2)&gapAUDphasesham<phasebinGap(gapBin2+1)))';
                    %  length(currGapbinATs)
                    dataGroupPool(subj,Session).hrbyStimbytACSbinsham(gapBin2,gapBin) = mean(trCurrGapTracksham(currGapbinATs));
                    dataGroupPool(subj,Session).hrbyStimbytACSbinshamNr(gapBin2,gapBin) = length(currGapbinATs);
                    
                end
            end
            for l = 1:length(lagBins)-1
                try
                    %fit by gap phase and stim_tacs phase lag
                    [dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,:), dataGroupPool(subj,Session).hrbytACSbinStimLagMain_RESNORM(l,:),dataGroupPool(subj,Session).hrbytACSbinStimLagMain_RESIDUAL(l,:),dataGroupPool(subj,Session).hrbytACSbinStimLagMain_prefPhase(l,:)] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbytACSbinStimLagMain(l,:) dataGroupPool(subj,Session).hrbytACSbinStimLagMain(l,:)]',X,params);
                    dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT_yAll(l,:) = dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,2) + dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,3).*(cos(X + dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,1))); % if you want to plot the predicted function
                    
                    [dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,:), dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESNORM(l,:),dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESIDUAL(l,:),dataGroupPool(subj,Session).hrbyStimbintACSLagMain_prefPhase(l,:)] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbyStimbintACSLagMain(l,:) dataGroupPool(subj,Session).hrbyStimbintACSLagMain(l,:)]',X,params);
                    dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT_yAll(l,:) = dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,2) + dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,1))); % if you want to plot the predicted function
                catch
                    errorPL = [subj Session l];
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(l,:) = nan;
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_RESNORM(l,:)= nan;
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_RESIDUAL(l,:)= nan;
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_prefPhase(l,:) = nan;
                dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT_yAll(l,:) = nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(l,:)= nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESNORM(l,:)= nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESIDUAL(l,:)= nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_prefPhase(l,:)= nan;
                dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT_yAll(l,:) = nan;

                end
            end
            %fit only by gap phas
            [dataGroupPool(subj,Session).hrbyStimbinMain_FIT, dataGroupPool(subj,Session).hrbyStimbinMain_RESNORM,dataGroupPool(subj,Session).hrbyStimbinMain_RESIDUAL,dataGroupPool(subj,Session).hrbyStimbinMain_prefPhase] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbyStimbinMain dataGroupPool(subj,Session).hrbyStimbinMain]',X,params);
            dataGroupPool(subj,Session).hrbyStimbinMain_FIT_yAll = dataGroupPool(subj,Session).hrbyStimbinMain_FIT(2) + dataGroupPool(subj,Session).hrbyStimbinMain_FIT(3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbinMain_FIT(1))); % if you want to plot the predicted function
            
            [dataGroupPool(subj,Session).hrbyStimbinsham_FIT, dataGroupPool(subj,Session).hrbyStimbinsham_RESNORM,dataGroupPool(subj,Session).hrbyStimbinsham_RESIDUAL,dataGroupPool(subj,Session).hrbyStimbinsham_prefPhase] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbyStimbinsham dataGroupPool(subj,Session).hrbyStimbinsham]',X,params);
            dataGroupPool(subj,Session).hrbyStimbinsham_FIT_yAll = dataGroupPool(subj,Session).hrbyStimbinsham_FIT(2) + dataGroupPool(subj,Session).hrbyStimbinsham_FIT(3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbinsham_FIT(1))); % if you want to plot the predicted function
            
            [dataGroupPool(subj,Session).hrbytACSbinMain_FIT, dataGroupPool(subj,Session).hrbytACSbinMain_RESNORM,dataGroupPool(subj,Session).hrbytACSbinMain_RESIDUAL,dataGroupPool(subj,Session).hrbytACSbinMain_prefPhase] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbytACSbinMain dataGroupPool(subj,Session).hrbytACSbinMain]',X,params);
            dataGroupPool(subj,Session).hrbytACSbinMain_FIT_yAll = dataGroupPool(subj,Session).hrbytACSbinMain_FIT(2) + dataGroupPool(subj,Session).hrbytACSbinMain_FIT(3).*(cos(X + dataGroupPool(subj,Session).hrbytACSbinMain_FIT(1))); % if you want to plot the predicted function
            
%             [dataGroupPool(subj,Session).hrbytACSbinsham_FIT, dataGroupPool(subj,Session).hrbytACSbinsham_RESNORM,dataGroupPool(subj,Session).hrbytACSbinsham_RESIDUAL,dataGroupPool(subj,Session).hrbytACSbinsham_prefPhase] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbytACSbinsham dataGroupPool(subj,Session).hrbytACSbinsham]',X,params);
%             dataGroupPool(subj,Session).hrbytACSbinsham_FIT_yAll = dataGroupPool(subj,Session).hrbytACSbinsham_FIT(2) + dataGroupPool(subj,Session).hrbytACSbinsham_FIT(3).*(cos(X + dataGroupPool(subj,Session).hrbytACSbinsham_FIT(1))); % if you want to plot the predicted function
            
            dataGroupPool(subj,Session).HRall = [dataGroupPool(subj,Session).hrbyStimbinsham_FIT(2) dataGroupPool(subj,Session).hrbyStimbinMain_FIT(2) dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(:,2)'];
%            dataGroupPool(subj,Session).HRalltacs = [dataGroupPool(subj,Session).hrbytACSbinsham_FIT(2) dataGroupPool(subj,Session).hrbytACSbinMain_FIT(:,2) dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(:,2)'];
            dataGroupPool(subj,Session).HRalltacs = [nan dataGroupPool(subj,Session).hrbytACSbinMain_FIT(2) dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(:,2)'];
            
            dataGroupPool(subj,Session).Ampall = [dataGroupPool(subj,Session).hrbyStimbinsham_FIT(3) dataGroupPool(subj,Session).hrbyStimbinMain_FIT(3) dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(:,3)'];
%            dataGroupPool(subj,Session).Ampalltacs = [dataGroupPool(subj,Session).hrbytACSbinsham_FIT(3) dataGroupPool(subj,Session).hrbytACSbinMain_FIT(:,3) dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(:,3)'];
            dataGroupPool(subj,Session).Ampalltacs = [nan dataGroupPool(subj,Session).hrbytACSbinMain_FIT(3) dataGroupPool(subj,Session).hrbytACSbinStimLagMain_FIT(:,3)'];
 
close all
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %create surrogate datasets by shuffeling the tACS phase lag
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gapAUDtACSlagMainS = nan(length(gapAUDtACSlagMain),NrIte);
            %               for ite =1:NrIte
            ite=0;
            while ite<NrIte+1
                try
                    ite=ite+1;
                    disp(ite)
                    neworder = randperm(length(gapAUDtACSlagMain));
                    gapAUDtACSlagMainS(:,ite) = gapAUDtACSlagMain(neworder);
                    %gapAUDtACSlagMain = wrapTo2Pi(gapAUD_tACSlag(gapStimcond==1)); %phase lag between tACS and audio signal
                    for gapBin = 1:length(phasebinGap)-1
                        currGapbinStim = find(gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))'; %gaps in current AUD bin
                        length(currGapbinStim)
                        for b=1:length(lagBins)-1 %separate by tACS-AUD phase lag
                            if lagBins1(b)>lagBins2(b)
                                currGapStimbinSlag = find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                                    (gapAUDtACSlagMainS(:,ite)>=lagBins1(b) | gapAUDtACSlagMainS(:,ite)<lagBins2(b)))';
                            else
                                currGapStimbinSlag =find((gapAUDphaseMain>=phasebinGap(gapBin)&gapAUDphaseMain<phasebinGap(gapBin+1))&...
                                    (gapAUDtACSlagMainS(:,ite)>=lagBins1(b)& gapAUDtACSlagMainS(:,ite)<lagBins2(b)))';
                            end
                            length(currGapStimbinSlag)
                            dataGroupPool(subj,Session).hrbyStimbintACSLagMainS(ite,b,gapBin) = mean(trCurrGapTrackMain(currGapStimbinSlag));
                            
%                                                     figure,
%                                                     subplot(1,2,1), circ_plot(gapAUDphaseMain(currGapStimbinSlag))
%                                                     title(['AUDIO gap bin: ' num2str(gapBin) ' lag bin' num2str(b)])
%                                                     subplot(1,2,2), circ_plot(gapAUDtACSlagMainS(currGapStimbinSlag,ite))
%                                                     title('lag')
                        end
                    end
                    for l = 1:length(lagBins)-1
                        %                       try
                        %fit by gap phase and stim_tacs phase lag
                        [dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,l,:), dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESNORMS(ite,l,:),dataGroupPool(subj,Session).hrbyStimbintACSLagMain_RESIDUALS(ite,l,:),dataGroupPool(subj,Session).hrbyStimbintACSLagMain_prefPhaseS(ite,l,:)] = fitCos2RatebyPhase ([squeeze(dataGroupPool(subj,Session).hrbyStimbintACSLagMainS(ite,l,:)); squeeze(dataGroupPool(subj,Session).hrbyStimbintACSLagMainS(ite,l,:))],X,params);
                        dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT_yAllS(ite,l,:) = dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,l,2) + dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,l,3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,l,1))); % if you want to plot the predicted function
                        %                       catch
                        %                       end
                    end
                    
                    %shuffle hit rates now and fit only by audio phase
                    trCurrGapTrackshamS = trCurrGapTracksham(randperm(length(trCurrGapTracksham)));
                    for gapBin = 1:length(phasebinGap)-1
                        currGapbinStimshamS = find(gapAUDphasesham>=phasebinGap(gapBin)&gapAUDphasesham<phasebinGap(gapBin+1))';
                        %    length(currGapbinStimshamS)
                        dataGroupPool(subj,Session).hrbyStimbinshamS(ite,gapBin) = mean(trCurrGapTrackshamS(currGapbinStimshamS));
                    end
                    %fit only by gap phas
                    [dataGroupPool(subj,Session).hrbyStimbinsham_FITS(ite,:), dataGroupPool(subj,Session).hrbyStimbinsham_RESNORMS(ite,:),dataGroupPool(subj,Session).hrbyStimbinsham_RESIDUALS(ite,:),dataGroupPool(subj,Session).hrbyStimbinsham_prefPhaseS(ite,:)] = fitCos2RatebyPhase ([dataGroupPool(subj,Session).hrbyStimbinshamS(ite,:) dataGroupPool(subj,Session).hrbyStimbinshamS(ite,:)]',X,params);
                    dataGroupPool(subj,Session).hrbyStimbinsham_FIT_yAllS(ite,:) = dataGroupPool(subj,Session).hrbyStimbinsham_FITS(ite,2) + dataGroupPool(subj,Session).hrbyStimbinsham_FITS(ite,3).*(cos(X + dataGroupPool(subj,Session).hrbyStimbinsham_FITS(ite,1))); % if you want to plot the predicted function
                    close all
                catch
                    ite =ite-1;
                end
            end
            %               %new strucut with multiple bins and overlaps
            %               dataLast(subj,Session).data = dataGroupPool(subj,Session);
        end
    end
end

save (fullfile(OUTdir,'dataGroupPooltACS_FMmodDepth_2024.mat'), 'dataGroupPool', '-v7.3')
load(fullfile(OUTdir,'dataGroupPooltACS_FMmodDepth_2024.mat'),'dataGroupPool')


%fit cos to tacs lag for session 2 as well
missingData = zeros(size(dataGroupPool));
fmlag = [];
fmlagFit = [];
for Session=1:2
    figure(Session)
    sp=1;
    for subj = 1:size(dataGroupPool,1)
       try
            groupmatrixFin.data(subj,:,Session) = dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FIT(:,3)';
            X         = 0:2*pi/NrAudtACSlag:2*pi-2*pi/NrAudtACSlag; %this is my phase vector
            Xplot     = (0:2*pi/1000:2*pi-2*pi/NrAudtACSlag)+2*pi/NrAudtACSlag/2; %this is my phase vector
            X         = X+2*pi/NrAudtACSlag/2;
            lag       = pi;
            intercept = 0.25;
            amp       = 0.1;
            params    = [lag intercept amp];
            
            TEMPdata = groupmatrixFin.data(subj,:,Session);
            [groupmatrixFin.dataFIT(subj,:,Session), ~,~,groupmatrixFin.dataprefPhase(subj,Session)] = fitCos2RatebyPhase ([TEMPdata TEMPdata]',X,params);
            groupmatrixFin.AlignCos_FIT_yAll(subj,Session,:) = groupmatrixFin.dataFIT(subj,2,Session) + groupmatrixFin.dataFIT(subj,3,Session).*(cos(Xplot + groupmatrixFin.dataFIT(subj,1,Session)))'; % if you want to plot the predicted function
            
            
         %   run the fit for the permuted datasets
            for ite =1:1000
                groupmatrixFin.dataS(subj,ite,:,Session)=dataGroupPool(subj,Session).hrbyStimbintACSLagMain_FITS(ite,:,3);
                TEMPdataS = squeeze(groupmatrixFin.dataS(subj,ite,:,Session))';
                [groupmatrixFin.dataFITS(subj,ite,:,Session), ~,~,groupmatrixFin.dataprefPhaseS(subj,ite,Session)] = fitCos2RatebyPhase ([TEMPdataS TEMPdataS]',X,params);
                groupmatrixFin.AlignCos_FIT_yAllS(subj,ite,:,Session) = groupmatrixFin.dataFITS(subj,ite,2,Session) + groupmatrixFin.dataFITS(subj,ite,3,Session).*(cos(X + groupmatrixFin.dataFITS(subj,ite,1,Session))); % if you want to plot the predicted function
            end
            subplot(7,5,subj)
            plot(X,squeeze(groupmatrixFin.data(subj,:,Session)))
            hold on
            plot(Xplot,squeeze(groupmatrixFin.AlignCos_FIT_yAll(subj,Session,:)))
            
            fmlag(subj,:,Session) = squeeze(groupmatrixFin.data(subj,:,Session));
            fmlagFit(subj,Session,:) = squeeze(groupmatrixFin.AlignCos_FIT_yAll(subj,Session,:));
            
            %ylim([0 0.4])
            sp=sp+1;
       catch
            missingData(subj,Session) = 1;
            groupmatrixFin.dataFIT(subj,:,Session)=nan;
            groupmatrixFin.dataprefPhase(subj,Session)=nan;
        end
        
    end
end
% filename = '/mnt/beegfs/users/yuranny.cabral/papertACSFMdepth/toShare/S8_Data_FigS2_Exp1.xlsx';
% writeXLS_FIgS2(filename,fmlag,fmlagFit);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%QUESTION 1: FM stimulus modulates behavior in sham?
Data2Plot = dataGroupPool;
Xnew = (0:2*pi/1000:2*pi-2*pi/NrAudtACSlag)+2*pi/NrAudtACSlag/2;
Xorig = 0:2*pi/9:2*pi-2*pi/9;
colorl = [247 148 189; 145 130 217; 100 158 166; 209 222 89; 247 195 74; 232 139 93; 27 118 188]./255;
origData = [];
origDataFit = [];
origDataSham = [];
origDataFitSham = [];

%%%FIGURE 1: SINGLE SUBJECT EXAMPLES
sp=1;    
for session = 1:2
    figure(session+100)
    for subj=1:length(Data2Plot)       
        %subplot(2,3,sp)
        try
        %(session)
        subplot(5,7,subj)
        for l =1:6
           % plot([Xorig Xorig+2*pi],[Data2Plot(subj,session).hrbyStimbintACSLagMain(l,:) Data2Plot(subj,session).hrbyStimbintACSLagMain(l,:)]','color',colorl(l,:))
            
            origData(session,subj,l,:) = Data2Plot(subj,session).hrbyStimbintACSLagMain(l,:);
            plot([Xorig Xorig+2*pi],[squeeze(origData(session,subj,l,:));squeeze(origData(session,subj,l,:))]','color',colorl(l,:))

            hold on
            fit = Data2Plot(subj,session).hrbyStimbintACSLagMain_FIT(l,:);
            y = fit(2) + fit(3).*(cos(Xnew + fit(1))); % if you want to plot the predicted function
            %plot([Xnew Xnew+2*pi],[y y],'--','color',colorl(l,:))
            origDataFit(session,subj,l,:) = y;
            plot([Xnew Xnew+2*pi],[squeeze(origDataFit(session,subj,l,:));squeeze(origDataFit(session,subj,l,:))],'--','color',colorl(l,:))
          %  scatter(Data2Plot(subj,session).hrbyStimbintACSLagMain_prefPhase(l),max(Data2Plot(subj,session).hrbyStimbintACSLagMain(l,:)))
        end
        origDataSham(session,subj,:) = Data2Plot(subj,session).hrbyStimbinsham ;
     %   plot([Xorig Xorig+2*pi],[Data2Plot(subj,session).hrbyStimbinsham Data2Plot(subj,session).hrbyStimbinsham],'color',colorl(7,:))
        plot([Xorig Xorig+2*pi],[ squeeze(origDataSham(session,subj,:));squeeze(origDataSham(session,subj,:))],'color',colorl(7,:))
        fit = Data2Plot(subj,session).hrbyStimbinsham_FIT;
        y = fit(2) + fit(3).*(cos(Xnew + fit(1))); % if you want to plot the predicted function
%         plot([Xnew Xnew+2*pi],[y y],'--','color',colorl(7,:))
        origDataFitSham(session,subj,:) = y;
        plot([Xnew Xnew+2*pi],[squeeze(origDataFitSham(session,subj,:)); squeeze(origDataFitSham(session,subj,:))],'--','color',colorl(7,:))
     %   scatter(Data2Plot(subj,session).hrbyStimbinsham_prefPhase,max(Data2Plot(subj,session).hrbyStimbinsham))
        title(['Subj_' num2str(subj) ' ModDepth ' num2str(session)])
        sp=sp+1;
                 catch
                 end
        ylim([0 1.2])
    end
    legend
end

% %convert variable to table and write to .xls file
% filename = '/mnt/beegfs/users/yuranny.cabral/papertACSFMdepth/toShare/S7_Data_FigS1_Exp1.xlsx';
% writeXLS_FIgS1(filename,origData, origDataSham,origDataFit, origDataFitSham)
% 
ampGroup           = nan(length(Data2Plot),2,7);
PhaseGroup         = nan(length(Data2Plot),2,7);
sigGroupSham       = nan(length(Data2Plot),2);
sigGroupShamPer95  = nan(length(Data2Plot),2);
sigGroupPerc95AmpS = nan(length(Data2Plot),2);
sigGroupPerc50AmpS = nan(length(Data2Plot),2);
meanAmpSbySession  = nan(length(Data2Plot),2);

for subj=1:length(Data2Plot)
    for session =1:2
        try
            ampGroup(subj,session,:) = [Data2Plot(subj,session).hrbyStimbinsham_FIT(3);Data2Plot(subj,session).hrbyStimbintACSLagMain_FIT(:,3)];
            HRGroup(subj,session,:) = [Data2Plot(subj,session).hrbyStimbinsham_FIT(2);Data2Plot(subj,session).hrbyStimbintACSLagMain_FIT(:,2)];

            PhaseGroup(subj,session,:) = [Data2Plot(subj,session).hrbyStimbinsham_prefPhase;Data2Plot(subj,session).hrbyStimbintACSLagMain_prefPhase];
            %estimate significance relative to perm distribution
            % for l=1:6
            % sigGroup(subj,session,l) = mean(Data2Plot(subj,session).data.hrbyStimbintACSLagMain_FIT(l,3)<squeeze(Data2Plot(subj,session).data.hrbyStimbintACSLagMain_FITS(:,l,3)));
            %
            % end
            sigGroupSham(subj,session)       = mean(Data2Plot(subj,session).hrbyStimbinsham_FIT(3)<squeeze(Data2Plot(subj,session).hrbyStimbinsham_FITS(:,3)));
            sigGroupShamPer95(subj,session)  = prctile(squeeze(Data2Plot(subj,session).hrbyStimbinsham_FITS(:,3)),95);
            sigGroupPerc95AmpS(subj,session) = prctile(squeeze(dataAmpS(subj,:,session)),95);
            sigGroupPerc50AmpS(subj,session) = prctile(squeeze(dataAmpS(subj,:,session)),50);
            meanAmpSbySession(subj,session)  = mean(squeeze(dataAmpS(subj,:,session)));
        catch
        end
    end
end

%SUBJECT TO EXCLUDE BECAUSE OF MISSING SESSION OR BINS
excSubj               = [1 6 8 9 12 21 26];
missingData           = ones(length(SUBJlist),1);
missingData (excSubj) = 0;
missingData           = logical(missingData);

ampGroupsameN           = ampGroup(missingData,:,:);
hrGroupsameN            = HRGroup(missingData,:,:);
PhaseGroupsameN         = PhaseGroup(missingData,:,:);
sigGroupShamsameN       = sigGroupSham(missingData,:);
sigGroupShamPerc95sameN = sigGroupShamPer95(missingData,:);

statsAll.nrSinMod2FM_sham = sum(sigGroupShamsameN<0.05);
sum(sigGroupShamsameN<9.2593e-04);%bonferroni correction
lillietest(ampGroupsameN(:,1))
lillietest(ampGroupsameN(:,2))

[statsAll.FMamp11vs39_H, statsAll.FMamp11vs39_pval, statsAll.FMamp11vs39_CI, statsAll.FMamp11vs39.stats] = ttest(ampGroupsameN(:,2,1),ampGroupsameN(:,1,1));
%statsAll.FMamp11vs39_cohenD = meanEffectSize(ampGroupsameN(:,2,1),ampGroupsameN(:,1,1),'cohen');

figure
distributionPlotYCC([ampGroupsameN(:,1,1),ampGroupsameN(:,2,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('Behavioral entrainment by modulation depth')
xlabel('modulation depth (%)')
xticklabels({'11' '39'})
ylabel('entAmp-FM')
ylim([-.05 .36])

% figure,scatter(ampGroupsameN(:,1,1),ampGroupsameN(:,2,1),'.')
% hold on
% plot([0 0.5],[0 0.5])
% title('Amplitude Sham')
% xlabel('11%')
% ylabel('39%')
% xlim([0 0.5])
% ylim([0 0.5])
% axis square

% %preferred phase
% for session =1:2
%     figure
%     for subj=1:length(PhaseGroupsameN)
%         subplot(4,8,subj)
%         circ_plot(squeeze(PhaseGroupsameN(subj,session,:)),'pretty')
%     end
% end

% figure,
% subplot(2,2,1)
% distributionPlot(ampGroupsameN(:,:,1), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2,'xyOri','flipped')
% hold on
% plot([mean(sigGroupShamPerc95sameN(:,1)) mean(sigGroupShamPerc95sameN(:,1))],[0 3])
% plot([mean(sigGroupShamPerc95sameN(:,2)) mean(sigGroupShamPerc95sameN(:,2))],[0 3],'--')

% xlim([0 0.5])
% title('FM-driven modulation')
% ylabel('ModDepth')
% xlabel('Amplitude')
% 
% subplot(2,2,2)
% circ_plot(PhaseGroupsameN(:,1,1),'pretty')
% title('S1')
% subplot(2,2,3)
% circ_plot(PhaseGroupsameN(:,2,1),'pretty')
% title('S2')
% subplot(2,2,4)
% circ_plot(circ_dist(PhaseGroupsameN(:,1,1),PhaseGroupsameN(:,2,1)),'pretty')
% title('Inter-ModDepth distance')

%Q2: is that reliable?
%compare amplitude and phase
[statsAll.S1vsS2ttestamp_H, statsAll.S1vsS2ttestamp_p, statsAll.S1vsS2ttestamp_CI, statsAll.S1vsS2ttestamp_stats]=ttest(ampGroupsameN(:,1,1),ampGroupsameN(:,2,1));
[statsAll.phaseclustShamP(1), statsAll.phaseclustShamZ(1)] = circ_rtest(PhaseGroupsameN(:,1,1));
[statsAll.phaseclustShamP(2), statsAll.phaseclustShamZ(2)] = circ_rtest(PhaseGroupsameN(:,2,1));
[statsAll.phaseclustSessionDistP, statsAll.phaseclustSessionDistV] = circ_vtest(circ_dist(PhaseGroupsameN(:,1,1),PhaseGroupsameN(:,2,1)),0);
[statsAll.phasecorrRho, statsAll.phasecorrP] = circ_corrcc(PhaseGroupsameN(:,1,1),PhaseGroupsameN(:,2,1));
[statsAll.ampcorrRho, statsAll.ampcorrP] = corr(ampGroupsameN(:,1,1),ampGroupsameN(:,2,1),'rows','complete');


% %Q3: does tacs modulates optimal stim phase?
% figure
% for l=1:6
%     subplot(2,6,l)
%     curPhase1 = PhaseGroupsameN(:,1,l+1);
%     [phase11p(l), phase11V(l)] = circ_vtest(curPhase1(~isnan(curPhase1)),circ_mean(PhaseGroupsameN(:,1,1)));
%     circ_plot(curPhase1(~isnan(curPhase1)),'pretty')
%     title('ModDepth 11%')
%     subplot(2,6,l+6)
%     curPhase2 = PhaseGroupsameN(:,2,l+1);
%     [phase39p(l), phase39V(l)] = circ_vtest(curPhase2(~isnan(curPhase2)),circ_mean(PhaseGroupsameN(:,2,1)));
%     circ_plot(curPhase2(~isnan(curPhase2)),'pretty')
%     title('ModDepth 39%')
% end
% writeXLS_FIg1(filename,origData, origDataSham,ampGroupsameN,gapSizeGroup_FMdepth,missingData,hrGroupsameN,fmlag,fmlagFit,NrFA)


%check if each phase lag condiftion is clustered around the same mean
for l=1:7
    for session=1:2
        [statsAll.phaseclustlagsP(l,session), statsAll.phaseclustlagsV(l,session)] = circ_vtest(PhaseGroupsameN(:,session,l),circ_mean(PhaseGroupsameN,[],'all'));
    end
end

%4. CHECK FOR DIFFERENCES IN GAP SIZE
%correctly move gap size for each mod depth instead of session
for s=1:length(SUBJlist)
    if SUBJlist{s,2}==11
        gapSizeGroup_FMdepth(s,:) = gapSizeGroup(s,:); 
    else %invert position
       gapSizeGroup_FMdepth(s,:) = [gapSizeGroup(s,2) gapSizeGroup(s,1)]; 
    end
end


[statsAll.S1vsS2ttestGAPsize_H, statsAll.S1vsS2ttestGAPsize_p, statsAll.S1vsS2ttestGAPsize_CI, statsAll.S1vsS2ttestGAPsize_stats]=ttest(gapSizeGroup_FMdepth(missingData,1),gapSizeGroup_FMdepth(missingData,2));
% lillietest(gapSizeGroup_FMdepth(:,1))
% lillietest(gapSizeGroup_FMdepth(:,2))
mean(gapSizeGroup_FMdepth(missingData,:),'all')
min(gapSizeGroup_FMdepth(missingData,:))
max(gapSizeGroup_FMdepth(missingData,:))
std(gapSizeGroup_FMdepth(missingData,:),[],'all')

figure
distributionPlot(gapSizeGroup_FMdepth(missingData,:),'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
title('GAP size')
xlabel('Gap size (ms)')
ylabel('ModDepth')

figure
distributionPlot([hrGroupsameN(:,1,1),hrGroupsameN(:,2,1)],'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
title('Hit rate')
xlabel('ModDepth')
ylabel('Hit rate')


%diff in HR
[statsAll.FMhr11vs39_H, statsAll.FMhr11vs39_pval, statsAll.FMhr11vs39_CI, statsAll.FMhr11vs39.stats] = ttest(hrGroupsameN(:,2,1),hrGroupsameN(:,1,1));
%statsAll.FMhr11vs39_cohenD = meanEffectSize(hrGroupsameN(:,2,1),hrGroupsameN(:,1,1),'cohen');

%Q5: does tACS modulates entrainment amplitude?

%align data using peak or optimal phase from cosine fit
ampReAligSubjCos1   = nan(length(SUBJlist),NrAudtACSlag+1,2);
ampReAligSubjRaw    = nan(length(SUBJlist),NrAudtACSlag+1,2);
ampReAligSubjCos1S   = nan(length(SUBJlist),1000,NrAudtACSlag+1,2);
ampReAligSubjRawS    = nan(length(SUBJlist),1000,NrAudtACSlag+1,2);

for session =1:2
    for subj=1:length(groupmatrixFin.dataprefPhase)
        if ~isnan(groupmatrixFin.data(subj,:,session))
            [optphaselag, optphaselagPOS] = min(abs(circ_dist(lagBins(1:6),groupmatrixFin.dataprefPhase(subj,session))));
            [optphaselagRaw, optphaselagPOSRaw] = max(groupmatrixFin.data(subj,:,session));
            ampReAligSubjCos1(subj,:,session) = [ampGroup(subj,session,1);wshift('1D',groupmatrixFin.data(subj,:,session),optphaselagPOS-1)'];
            ampReAligSubjRaw(subj,:,session) = [ampGroup(subj,session,1);wshift('1D',groupmatrixFin.data(subj,:,session),optphaselagPOSRaw-1)'];
        end
        %do the same for the surrogate data
        for ite =1:1000
             if ~isnan(groupmatrixFin.dataS(subj,ite,:,session))
            [optphaselagS, optphaselagPOSS] = min(abs(circ_dist(lagBins(1:6),groupmatrixFin.dataprefPhaseS(subj,ite,session))));
            [optphaselagRawS, optphaselagPOSRawS] = max(groupmatrixFin.dataS(subj,ite,:,session));
            ampReAligSubjCos1S(subj,ite,:,session) = [ampGroup(subj,session,1);wshift('1D',squeeze(groupmatrixFin.dataS(subj,ite,:,session))',optphaselagPOSS-1)'];
            ampReAligSubjRawS(subj,ite,:,session) = [ampGroup(subj,session,1);wshift('1D',squeeze(groupmatrixFin.dataS(subj,ite,:,session))',optphaselagPOSRawS-1)'];
             end
        end
    end
end
figure, plot(ampReAligSubjCos1(subj,2:7,2))
hold on
plot (groupmatrixFin.data(subj,:,2))
plot(optphaselagPOS, max(groupmatrixFin.AlignCos_FIT_yAll(subj,2,:)),'*')
plot(Xplot,squeeze(groupmatrixFin.AlignCos_FIT_yAll(subj,2,:)))

% figure
% subplot(2,2,1)
% distributionPlotYCC(squeeze(ampReAligSubjCos1(missingData,:,1)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
% 
% subplot(2,2,2)
% distributionPlotYCC(squeeze(ampReAligSubjCos1(missingData,:,2)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
% 
% subplot(2,2,3)
% distributionPlotYCC(squeeze(ampReAligSubjRaw(missingData,:,1)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
% 
% subplot(2,2,4)
% distributionPlotYCC(squeeze(ampReAligSubjRaw(missingData,:,2)), 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')


maxadj = squeeze(mean(ampReAligSubjCos1(:,[3 7],:),2));
minadj = squeeze(mean(ampReAligSubjCos1(:,[4 6],:),2));

maxadjS = squeeze(mean(ampReAligSubjCos1S(:,:,[3 7],:),3));
minadjS = squeeze(mean(ampReAligSubjCos1S(:,:,[4 6],:),3));
   
%get subjMod diff from sham for real and surrogate
maxadj2Sham = maxadj-squeeze(ampReAligSubjCos1(:,1,:));
minadj2Sham = minadj-squeeze(ampReAligSubjCos1(:,1,:));
maxadj2min = maxadj-minadj;

maxadj2ShamS = maxadjS-squeeze(ampReAligSubjCos1S(:,:,1,:));
minadj2ShamS = minadjS-squeeze(ampReAligSubjCos1S(:,:,1,:));
maxadj2minS = maxadjS-minadjS;

%compute z-scores relative to surrogate distribution
maxadj2Sham_zscore2S = (maxadj2Sham-squeeze(mean(maxadj2ShamS,2)))./squeeze(std(maxadj2ShamS,[],2));
minadj2Sham_zscore2S = (minadj2Sham-squeeze(mean(minadj2ShamS,2)))./squeeze(std(minadj2ShamS,[],2));
maxadj2min_zscore2S = (maxadj2min-squeeze(mean(maxadj2minS,2)))./squeeze(std(maxadj2minS,[],2));

statsAll.MaxAdjanova_rm_S1 = anova_rm([squeeze(ampReAligSubjCos1(missingData,1,1)) maxadj(missingData,1) minadj(missingData,1)]);
statsAll.MaxAdjanova_rm_S2 = anova_rm([squeeze(ampReAligSubjCos1(missingData,1,2)) maxadj(missingData,2) minadj(missingData,2)]);
[statsAll.MaxAdjanova_rm_All,statsAll.MaxAdjanova_rm_AllTable] = anova_rm({[squeeze(ampReAligSubjCos1(missingData,1,1)) maxadj(missingData,1) minadj(missingData,1)] [squeeze(ampReAligSubjCos1(missingData,1,2)) maxadj(missingData,2) minadj(missingData,2)]});
[statsAll.posthocttest_S1_H(1), statsAll.posthocttest_S1_P(1), statsAll.posthocttest_S1_CI(1,:), statsAll.posthocttest_S1_Stats(1).stats]=ttest(squeeze(ampReAligSubjCos1(missingData,1,1)),maxadj(missingData,1));
[statsAll.posthocttest_S1_H(2), statsAll.posthocttest_S1_P(2), statsAll.posthocttest_S1_CI(2,:), statsAll.posthocttest_S1_Stats(2).stats]=ttest(squeeze(ampReAligSubjCos1(missingData,1,1)),minadj(missingData,1));
[statsAll.posthocttest_S1_H(3), statsAll.posthocttest_S1_P(3), statsAll.posthocttest_S1_CI(3,:), statsAll.posthocttest_S1_Stats(3).stats]=ttest(maxadj(missingData,1),minadj(missingData,1));
[statsAll.posthocttest_S2_H(1), statsAll.posthocttest_S2_P(1), statsAll.posthocttest_S2_CI(1,:), statsAll.posthocttest_S2_Stats(1).stats]=ttest(squeeze(ampReAligSubjCos1(missingData,1,2)),maxadj(missingData,2));
[statsAll.posthocttest_S2_H(2), statsAll.posthocttest_S2_P(2), statsAll.posthocttest_S2_CI(2,:), statsAll.posthocttest_S2_Stats(2).stats]=ttest(squeeze(ampReAligSubjCos1(missingData,1,2)),minadj(missingData,2));
[statsAll.posthocttest_S2_H(3), statsAll.posthocttest_S2_P(3), statsAll.posthocttest_S2_CI(3,:), statsAll.posthocttest_S2_Stats(3).stats]=ttest(maxadj(missingData,2),minadj(missingData,2));
statsAll.posthocttest_Bonferroni = [statsAll.posthocttest_S1_P*6 statsAll.posthocttest_S2_P*6];

tabletACS = table(squeeze(ampReAligSubjCos1(missingData,1,1)), maxadj(missingData,1), minadj(missingData,1), squeeze(ampReAligSubjCos1(missingData,1,2)), maxadj(missingData,2), minadj(missingData,2), (1:length(minadj(missingData,2)))','VariableNames',...
    {'sham1','pos1','neg1','sham2','pos2', 'neg2','subjID'});

withinDesign = table(([1 2 3 1 2 3])',([1 1 1 2 2 2])','VariableNames',{'tACS','ModDepth'});
withinDesign.tACS = categorical(withinDesign.tACS);
withinDesign.ModDepth = categorical(withinDesign.ModDepth);
%rmfitting - table - testing groups - grouping
myRM = fitrm(tabletACS,'sham1-neg2 ~ 1','withinDesign',withinDesign);
 
[b1_rm,A_rm,C_rm,D_rm] = ranova(myRM,'WithinModel','tACS*ModDepth');
Multicompare_tableTACS = multcompare(myRM,'tACS');
Multicompare_tableMD = multcompare(myRM,'ModDepth');

figure, plotprofile(myRM,'tACS')
figure, plotprofile(myRM,'ModDepth')

figure
subplot(1,2,1)
distributionPlotYCC([ampReAligSubjCos1(missingData,1,1) ampReAligSubjCos1(missingData,7,1) ampReAligSubjCos1(missingData,2:6,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('FM-driven modulation S1')
ylim([-0.1 0.6])
xlabel('Realigned tACS lags')
xticklabels({'sham' '-pi/3' '0' 'pi/3' '2pi/3' '+-pi' '-2pi/3'})
ylabel('Amplitude')
subplot(1,2,2)
distributionPlotYCC([ampReAligSubjCos1(missingData,1,2) ampReAligSubjCos1(missingData,7,2) ampReAligSubjCos1(missingData,2:6,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('FM-driven modulation S2')
ylim([-0.1 0.6])
xlabel('Realigned tACS lags')
ylabel('Amplitude')
xticklabels({'sham' '-pi/3' '0' 'pi/3' '2pi/3' '+-pi' '-2pi/3'})

figure
subplot(1,2,1)
distributionPlotYCC([squeeze(ampReAligSubjCos1(missingData,1,1)) maxadj(missingData,1) minadj(missingData,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-0.1 0.5])
title('Realigned amplitude 11 MD')
ylabel('Amplitude')
xlabel('Condition')
xticklabels({'sham' 'pos' 'neg'})
subplot(1,2,2)
distributionPlotYCC([squeeze(ampReAligSubjCos1(missingData,1,2)) maxadj(missingData,2) minadj(missingData,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-0.1 0.5])
title('Realigned amplitude 39 MD')
ylabel('Amplitude')
xlabel('Condition')
xticklabels({'sham' 'pos' 'neg'})

%PLOT z-score by Surrogate
%remove missing values
maxadj2Sham_zscore2S = maxadj2Sham_zscore2S(missingData,:);
minadj2Sham_zscore2S = minadj2Sham_zscore2S(missingData,:);
maxadj2min_zscore2S = maxadj2min_zscore2S(missingData,:);

figure
subplot(1,2,1)
distributionPlotYCC([maxadj2Sham_zscore2S(:,1) minadj2Sham_zscore2S(:,1) maxadj2min_zscore2S(:,1)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-3 3])
title('Normalized Realigned amplitude 11 MD')
ylabel('Amplitude (z-score)')
xlabel('Condition')
xticklabels({'max-sham' 'min-sham' 'max-min'})
subplot(1,2,2)
distributionPlotYCC([maxadj2Sham_zscore2S(:,2) minadj2Sham_zscore2S(:,2) maxadj2min_zscore2S(:,2)], 'showMM',1,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)
ylim([-3 3])
title('Normalized Realigned amplitude 39 MD')
ylabel('Amplitude (z-score)')
xlabel('Condition')
xticklabels({'max-sham' 'min-sham' 'max-min'})

%ttest against zero to check for effect different than chance
[statsAll.posthocttest_zscore_H(1), statsAll.posthocttest_zscore_P(1), statsAll.posthocttest_zscore_CI(1,:), statsAll.posthocttest_zscore_Stats(1).stats]=ttest(maxadj2Sham_zscore2S(:,1));
[statsAll.posthocttest_zscore_H(2), statsAll.posthocttest_zscore_P(2), statsAll.posthocttest_zscore_CI(2,:), statsAll.posthocttest_zscore_Stats(2).stats]=ttest(minadj2Sham_zscore2S(:,1));
[statsAll.posthocttest_zscore_H(3), statsAll.posthocttest_zscore_P(3), statsAll.posthocttest_zscore_CI(3,:), statsAll.posthocttest_zscore_Stats(3).stats]=ttest(maxadj2min_zscore2S(:,1));

[statsAll.posthocttest_zscore_H(4), statsAll.posthocttest_zscore_P(4), statsAll.posthocttest_zscore_CI(4,:), statsAll.posthocttest_zscore_Stats(4).stats]=ttest(maxadj2Sham_zscore2S(:,2));
[statsAll.posthocttest_zscore_H(5), statsAll.posthocttest_zscore_P(5), statsAll.posthocttest_zscore_CI(5,:), statsAll.posthocttest_zscore_Stats(5).stats]=ttest(minadj2Sham_zscore2S(:,2));
[statsAll.posthocttest_zscore_H(6), statsAll.posthocttest_zscore_P(6), statsAll.posthocttest_zscore_CI(6,:), statsAll.posthocttest_zscore_Stats(6).stats]=ttest(maxadj2min_zscore2S(:,2));

FMdepths = [11 39];
gapDetectALL = [];
gaptACSangleALL = [];
gapAUDangleALL = [];
subjID = [];
modDepth = [];

for subj =1:34
    if missingData(subj)==1
    for session = 1:2
        gapDetectALL      = [gapDetectALL; allDataModDepth(subj,Session).gapDetect(allDataModDepth(subj,Session).gapStimcond==1)];
        gaptACSangleALL   = [gaptACSangleALL; allDataModDepth(subj,Session).gaptACSangle(allDataModDepth(subj,Session).gapStimcond==1)'];
        gapAUDangleALL    = [gapAUDangleALL;allDataModDepth(subj,Session).gapAUDangle(allDataModDepth(subj,Session).gapStimcond==1)'];
        subjID            = [subjID; zeros(length(allDataModDepth(subj,Session).gapDetect(allDataModDepth(subj,Session).gapStimcond==1)),1)+subj];
        modDepth          = [modDepth; zeros(length(allDataModDepth(subj,Session).gapDetect(allDataModDepth(subj,Session).gapStimcond==1)),1)+FMdepths(session)];
    end
        
    end
end
missingVals = isnan(gapAUDangleALL);
gapDetectALL(missingVals)    = [];     
gaptACSangleALL(missingVals) = [];
gapAUDangleALL(missingVals)  = [];
subjID(missingVals)          = [];
modDepth(missingVals)        = [];
        
    mfr = table(gapDetectALL,zscore(sin(gaptACSangleALL)),zscore(cos(gaptACSangleALL)),zscore(sin(gapAUDangleALL)),zscore(cos(gapAUDangleALL)),subjID,modDepth,'VariableNames',...
    {'response','sintACS','costACS','sinAUD','cosAUD','subjectNr','modDepth'});

glme_1 = fitglme(mfr,'response ~ 1  + sinAUD + cosAUD + (1|subjectNr)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','effects','DispersionFlag',true);

glme_2 = fitglme(mfr,'response ~ 1  + sinAUD + cosAUD + modDepth + modDepth:sinAUD + modDepth:cosAUD +(1|subjectNr)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','effects','DispersionFlag',true);

glme_3 = fitglme(mfr,'response ~ 1  + sintACS + costACS +(1|subjectNr) + (sintACS|subjectNr) + (costACS|subjectNr)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','effects','DispersionFlag',true);

glme_4 = fitglme(mfr,'response ~ 1  + sinAUD + cosAUD + sintACS + costACS +(1|subjectNr) + (sintACS|subjectNr) + (costACS|subjectNr)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','effects','DispersionFlag',true);
glme_5 = fitglme(mfr,'response ~ 1  + sinAUD + cosAUD + sintACS + costACS + modDepth + modDepth:sintACS + modDepth:costACS + (1|subjectNr) + (sintACS|subjectNr) + (costACS|subjectNr)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','effects','DispersionFlag',true);

AICtable = [glme_1.ModelCriterion.AIC; glme_2.ModelCriterion.AIC;glme_3.ModelCriterion.AIC;glme_4.ModelCriterion.AIC;glme_5.ModelCriterion.AIC];
AICtable(:,2)=AICtable(:,1)-min(AICtable(:,1));

%WINNING MODEL = GLME1
%%run model for surrogate data to compare the auc
for ite=1:1000
    disp(ite)
    newDat = gapAUDangleALL;
    
    for subj = 1:max(subjID)
        currDat  = find(subjID==subj);
        if ~isempty(currDat)
        neworder = randperm(length(currDat))';
        tempDat  = gapAUDangleALL(currDat);
        newDat(currDat) = tempDat(neworder);
        end
    end
    disp(newDat)
    
    mfr_S = table(gapDetectALL,zscore(sin(gaptACSangleALL)),zscore(cos(gaptACSangleALL)),zscore(sin(newDat)),zscore(cos(newDat)),subjID,modDepth,'VariableNames',...
    {'response','sintACS','costACS','sinAUD','cosAUD','subjectNr','modDepth'});

glme_1S(1,ite).model = fitglme(mfr_S,'response ~ 1  + sinAUD + cosAUD + (1|subjectNr)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','effects','DispersionFlag',true);

scores = glme_1S(1,ite).model.Fitted;
[X1s,Y1s,T1s,AUC1s(ite,:)] = perfcurve(gapDetectALL,scores',1);
end

scores = glme_1.Fitted;
[X,Y,T,AUC] = perfcurve(gapDetectALL,scores',1);

AUC

figure
plot(X,Y)
hold on
plot(X1s,Y1s)
xlabel('False positive rate')
ylabel('True positive rate')
title('ROC for Classification by Logistic Regression Session 1')

%COMPUTE P-VALS FROM DIST OF AUC FROM SURROGATE DATASET
pval(1) = sum(AUC1s>AUC)./1000;

figure
histogram(AUC1s)
hold on
plot([AUC AUC],[0 180])
xlabel('AUC')
ylabel('# iterations')
title('AUC from simulated vs real data Session 1')
% ylim([0 180])
% xlim([0.664 0.676])

save tACSAna_MaxvsMin_FMmodDepth_April2024 -v7.3