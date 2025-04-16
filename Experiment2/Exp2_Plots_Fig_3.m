%FIGURES FOR PUBLICATION
load('Fig3_Exp2') 

%Figure 3b
%Effect of tACS frequency on HR and reaction time
nBins = 6;
X         = 0:2*pi/nBins:2*pi-2*pi/nBins; %this is my phase vector
X         = X+2*pi/nBins/2;
figure, subplot(1,2,1)
distributionPlot(GroupMeanHR,'color',[1 1 1], 'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(GroupMeanHR,'xValues',[1.3 2.3 3.3 4.3 5.3], 'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

title('mean HIT RATE per tACS condition')
xlabel('tACS condition (Hz)')
xticklabels({'Sham' '0.8' '2' '3.2' '4.4'})
ylabel('HIT RATE')
subplot(1,2,2)
distributionPlot(GroupMeanRT,'color',[1 1 1], 'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(GroupMeanRT,'xValues',[1.3 2.3 3.3 4.3 5.3], 'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('mean REACTION TIME per tACS condition')
xlabel('tACS condition (Hz)')
xticklabels({'Sham' '0.8' '2' '3.2' '4.4'})
ylabel('Reaction Time (s)')

%Figure 3d: Single subject data
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
    xlabel ('Phase')
    ylabel ('Hit rate')
    
    set(gca, 'XTick',X)
    xticklabels({'-pi' ,'', '', '', '', 'pi'})
    title (['HR by tACS phase P' num2str(s)])
end
legend({'0.8 Hz';''; '2';''; '3.2';''; '4.4'})
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',11,'FontWeight','Bold', 'LineWidth', 2);

%compute percent change

percentChange = [(tabletACS.tACS08-tabletACS.sham08)./tabletACS.sham08,...
    (tabletACS.tACS2-tabletACS.sham2)./tabletACS.sham2,...
    (tabletACS.tACS32-tabletACS.sham32)./tabletACS.sham32,...
    (tabletACS.tACS44-tabletACS.sham44)./tabletACS.sham44];

%Figure 3e   
figure
distributionPlot(percentChange,'color',[1 1 1], 'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
hold on
distributionPlot(percentChange,'xValues',[1.3 2.3 3.3 4.3], 'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('mod Amp by Freq')
xlabel('frequency')
set(gca, 'XTick',[0.8 2 3.2 4.4])
xticklabels({'0.8' ,'2', '3.2', '4.4'})
ylabel('amplitude')
