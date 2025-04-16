
load('Fig5_Exp3')
%Test for linear trend
ampTable = table([groupNormAmp(:,1); groupNormAmp(:,2); groupNormAmp(:,3); groupNormAmp(:,4)], [ones(length(groupNormAmp(:,1)),1)-.2; ones(length(groupNormAmp(:,1)),1)+1; ones(length(groupNormAmp(:,1)),1)+2.2; ones(length(groupNormAmp(:,1)),1)+3.4],'VariableNames',...
    {'Amp','Freq'});
Amplm = fitlm(ampTable,'Amp~Freq');
figure, plot(Amplm)
hold on
distributionPlot(groupNormAmp,'color',[1 1 1],'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
distributionPlot(groupNormAmp,'xValues',[1.3 2.3 3.3 4.3],'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
legend off
ylabel('Normalized amplitude')
xlabel('FM(Hz)')
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',11,'FontWeight','Bold', 'LineWidth', 2);

%VL
vlTable = table([groupNormVL(:,1); groupNormVL(:,2); groupNormVL(:,3); groupNormVL(:,4)], [ones(length(groupNormAmp(:,1)),1)-.2; ones(length(groupNormAmp(:,1)),1)+1; ones(length(groupNormAmp(:,1)),1)+2.2; ones(length(groupNormAmp(:,1)),1)+3.4],'VariableNames',...
    {'VL','Freq'});
VLlm = fitlm(vlTable,'VL~Freq');
figure, plot(VLlm)
hold on
distributionPlot(groupNormVL,'color',[1 1 1],'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
distributionPlot(groupNormVL,'xValues',[1.3 2.3 3.3 4.3],'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
legend off
ylabel('Normalized Vector length')
xlabel('FM(Hz)')
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',11,'FontWeight','Bold', 'LineWidth', 2);
