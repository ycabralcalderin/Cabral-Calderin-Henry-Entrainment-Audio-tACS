%Panel 1d
colorl = [247 148 189; 145 130 217; 100 158 166; 209 222 89; 247 195 74; 232 139 93; 27 118 188]./255;
t = 1/1000:1/1000:1-1/1000;
t2plot = (FMphase(end)+2*pi)/1000:(FMphase(end)+2*pi)/1000:(FMphase(end)+2*pi)-((FMphase(end)+2*pi)/1000);
f = 2;
ph = 0;
for session = 1:2
    h = figure(1);
    set(h, 'NumberTitle', 'off', ...
        'Name', 'Fig 1D');
    if session == 1
        for subj = 1:3
            subplot(2,3,subj)
            hold on
            plot(t2plot,normalize(cos(f*2*pi.*t),'range'),'k','linewidth',2)
            
            for l =1:6
                plot([FMphase FMphase+2*pi],[squeeze(HRbyFMphase_tACS_panelD(session,subj,l,:));squeeze(HRbyFMphase_tACS_panelD(session,subj,l,:))]','color',colorl(l,:),'linewidth',2)
            end
            plot([FMphase FMphase+2*pi],[squeeze(HRbyFMphase_Sham_panelD(session,subj,:));squeeze(HRbyFMphase_Sham_panelD(session,subj,:))]','color',colorl(7,:),'linewidth',2)
            
            ylim([0 1])
            yticks([0 0.5 1])
            yticklabels({'','',''})
            xticks([0 5 10 12])
            xticklabels({'','','',''})
        end
    else
        for subj = 1:3
            subplot(2,3,subj+3)
            plot(t2plot,normalize(cos(f*2*pi.*t),'range'),'k','linewidth',2)
            hold on
            
            for l =1:6
                plot([FMphase FMphase+2*pi],[squeeze(HRbyFMphase_tACS_panelD(session,subj,l,:));squeeze(HRbyFMphase_tACS_panelD(session,subj,l,:))]','color',colorl(l,:),'linewidth',2)
            end
            plot([FMphase FMphase+2*pi],[squeeze(HRbyFMphase_Sham_panelD(session,subj,:));squeeze(HRbyFMphase_Sham_panelD(session,subj,:))]','color',colorl(7,:),'linewidth',2)
            ylim([0 1])
            yticks([0 0.5 1])
            yticklabels({'','',''})
            xticks([0 5 10 12])
            xticklabels({'','','',''})
            hax = axes;
            yyaxis(hax,'right');
            hax.YAxis(2).Visible='off';
            
        end
    end
end
subplot(2,3,1)
ylabel('Hit rate','FontSize',12)
yticks([0 0.5 1])
yticklabels({'0','0.5','1'})

subplot(2,3,2)
title('11%', 'FontSize', 12)
subplot(2,3,4)
yticks([0 0.5 1])
yticklabels({'0','0.5','1'})
subplot(2,3,5)
title('39%', 'FontSize', 12)
xlabel('FM-phase','FontSize',12)
xticks([0 5 10 11.87])
xticklabels({'-2pi','','','+2pi'})
set(gcf, 'color', 'white');


%figure 1e
figure
distributionPlot(entAmp_FMSham,'color',[1 1 1],'xValues',[1 2], 'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(entAmp_FMSham,'xValues',[1.3 2.3], 'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('Behavioral entrainment by modulation depth')
xlabel('Modulation depth (%)')
xticks([1 2])
xticklabels({'11' '39'})
ylabel('entAmp-FM')
ylim([-.05 .36])

%figure 1f
figure
distributionPlot(gapSize,'color',[1 1 1],'xValues',[1 2], 'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(gapSize,'xValues',[1.3 2.3], 'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('Modulation depth (%)')
xticks([1 2])
xticklabels({'11' '39'})
ylabel('Gap size (ms)')
ylim([0 25])

%figure 1g
figure
distributionPlot(hrSham,'color',[1 1 1],'xValues',[1 2], 'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(hrSham,'xValues',[1.3 2.3], 'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('Modulation depth (%)')
xticks([1 2])
xticklabels({'11' '39'})
ylabel('Hit rate')
ylim([0 1])

%figure 1h
figure
distributionPlot(NrFA,'color',[1 1 1],'xValues',[1 2], 'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(NrFA,'xValues',[1.3 2.3], 'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('Modulation depth (%)')
xticks([1 2])
xticklabels({'11' '39'})
ylabel('FA')
%ylim([-10 65])

%panel 1i entAmp-FM by tACS lag
%single subject examples 
entAmp_FM_tACSlag = fmlag(subjs,:,:);
entAmp_FM_tACSlag_cosFit = fmlagFit(subjs,:,:);

tACSlags = X;
interp_tACSlags = Xplot;

figure
idx =1;

for subj = 1:3
    for session = 1:2
        h = figure(1);
        set(h, 'NumberTitle', 'off', ...
            'Name', 'Fig 1i'); 
        subplot(3,2,idx)
        hold on
        plot(tACSlags,squeeze(entAmp_FM_tACSlag(subj,:,session)),'k','linewidth',2)
        plot(interp_tACSlags,squeeze(entAmp_FM_tACSlag_cosFit(subj,session,:)),'--','color',[0.5 0.5 0.5],'linewidth',2)
        xticks([1 2 3 4 5 6])
        xticklabels({'0','pi/3','2pi/3','pi','-2pi/3','-pi/3'})
        idx=idx+1;
    end
    
end
subplot(3,2,1)
ylabel('entAmp-FM','FontSize',12)
title('11%', 'FontSize', 12)

subplot(3,2,2)
title('39%', 'FontSize', 12)

subplot(3,2,5)
xlabel('tACS lag','FontSize',12)
set(gcf, 'color', 'white');

%%--------------------
%Figure 2b
figure
subplot(1,2,1)
distributionPlot(entAmp_FM_bytACSlag_realigned_modDepth11,'color',[1 1 1],'xValues',[1 2 3 4 5 6 7], 'showMM',2,'addSpread',0,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(entAmp_FM_bytACSlag_realigned_modDepth11,'color',[1 1 1],'xValues',[0.5 1.5 2.5 3.5 4.5 5.5 6.5], 'showMM',0,'addSpread',1,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('Realigned tACS lags')
xticks([1 2 3 4 5 6 7])
xticklabels(tACSConds)
ylabel('entAmp-FM')
ylim([-.05 .4])
title('11%', 'FontSize',12)

subplot(1,2,2)
distributionPlot(entAmp_FM_bytACSlag_realigned_modDepth39,'color',[1 1 1],'xValues',[1 2 3 4 5 6 7], 'showMM',2,'addSpread',0,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(entAmp_FM_bytACSlag_realigned_modDepth39,'color',[1 1 1],'xValues',[0.5 1.5 2.5 3.5 4.5 5.5 6.5], 'showMM',0,'addSpread',1,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('Realigned tACS lags')
xticks([1 2 3 4 5 6 7])
xticklabels(tACSConds)
ylabel('entAmp-FM')
ylim([-.05 .4])
title('39%', 'FontSize',12)

%Figure 2d
figure
subplot(1,2,1)
distributionPlot(maxAdjdata_modDepth11,'color',[1 1 1],'xValues',[1 2 3], 'showMM',2,'addSpread',0,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(maxAdjdata_modDepth11,'color',[1 1 1],'xValues',[0.5 1.5 2.5], 'showMM',0,'addSpread',1,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('tACS condition')
xticks([1 2 3])
xticklabels(maxAdjdat_conditions)
ylabel('entAmp-FM')
ylim([-.05 .4])
title('11%', 'FontSize',12)

subplot(1,2,2)
distributionPlot(maxAdjdata_modDepth39,'color',[1 1 1],'xValues',[1 2 3], 'showMM',2,'addSpread',0,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(maxAdjdata_modDepth39,'color',[1 1 1],'xValues',[0.5 1.5 2.5], 'showMM',0,'addSpread',1,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('tACS condition')
xticks([1 2 3])
xticklabels(maxAdjdat_conditions)
ylabel('entAmp-FM')
ylim([-.05 .4])
title('39%', 'FontSize',12)

%Figure 2f
figure
subplot(1,2,1)
distributionPlot(maxAdjdata_modDepth11_zscore,'color',[1 1 1],'xValues',[1 2 3], 'showMM',2,'addSpread',0,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(maxAdjdata_modDepth11_zscore,'color',[1 1 1],'xValues',[0.5 1.5 2.5], 'showMM',0,'addSpread',1,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('tACS condition')
xticks([1 2 3])
xticklabels(maxAdjdat_conditions_zscore)
ylabel('entAmp-FM(z-score)')
ylim([-3 3])
title('11%', 'FontSize',12)

subplot(1,2,2)
distributionPlot(maxAdjdata_modDepth39_zscore,'color',[1 1 1],'xValues',[1 2 3], 'showMM',2,'addSpread',0,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(maxAdjdata_modDepth39_zscore,'color',[1 1 1],'xValues',[0.5 1.5 2.5], 'showMM',0,'addSpread',1,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('tACS condition')
xticks([1 2 3])
xticklabels(maxAdjdat_conditions_zscore)
ylabel('entAmp-FM(z-score)')
ylim([-3 3])
title('39%', 'FontSize',12)