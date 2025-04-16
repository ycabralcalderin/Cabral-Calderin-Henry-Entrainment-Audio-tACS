load('Fig6_Exp3')

%figure 6b
gapSizeTable = table([gapSize(:,1); gapSize(:,2); gapSize(:,3); gapSize(:,4)], [ones(length(gapSize(:,1)),1)-.2; ones(length(gapSize(:,1)),1)+1; ones(length(gapSize(:,1)),1)+2.2; ones(length(gapSize(:,1)),1)+3.4],'VariableNames',...
    {'gapSize','Freq'});
gapSizelm = fitlm(gapSizeTable,'gapSize~Freq');
figure, plot(gapSizelm,'linewidth',2)
hold on
distributionPlot(gapSize,'color',[1 1 1],'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
distributionPlot(gapSize,'xValues',[1.3 2.3 3.3 4.3],'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
legend off
ylabel('Gap size (ms)')
xlabel('FM(Hz)')
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',11,'FontWeight','Bold', 'LineWidth', 2);

%figure 6c
hRTable = table([hR(:,1); hR(:,2); hR(:,3); hR(:,4)], [ones(length(hR(:,1)),1)-.2; ones(length(hR(:,1)),1)+1; ones(length(hR(:,1)),1)+2.2; ones(length(hR(:,1)),1)+3.4],'VariableNames',...
    {'hR','Freq'});
hRlm = fitlm(hRTable,'hR~Freq');
figure, plot(hRlm)
hold on
distributionPlot(hR,'color',[1 1 1],'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
distributionPlot(hR,'xValues',[1.3 2.3 3.3 4.3],'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
legend off
ylabel('Hit rate')
xlabel('FM(Hz)')
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',11,'FontWeight','Bold', 'LineWidth', 2);

%figure 6d
groupfitampTable = table([groupfitamp(:,1); groupfitamp(:,2); groupfitamp(:,3); groupfitamp(:,4)], [ones(length(groupfitamp(:,1)),1)-.2; ones(length(groupfitamp(:,1)),1)+1; ones(length(groupfitamp(:,1)),1)+2.2; ones(length(groupfitamp(:,1)),1)+3.4],'VariableNames',...
    {'groupfitamp','Freq'});
groupfitamplm = fitlm(groupfitampTable,'groupfitamp~Freq');
figure, plot(groupfitamplm)
hold on
distributionPlot(groupfitamp,'color',[1 1 1],'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
distributionPlot(groupfitamp,'xValues',[1.3 2.3 3.3 4.3],'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
legend off
ylabel('EntAmp-FM')
xlabel('FM(Hz)')
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',11,'FontWeight','Bold', 'LineWidth', 2);

