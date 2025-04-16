load('Fig4_Exp2') 

%Figure 4d
figure 
distributionPlot(GroupRF','color',[1 1 1],'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
title('predicted Resonance Frequency')
ylabel('Predicted Resonance Frequency (Hz)','FontSize',12)
set(findobj(gcf,'type','axes'),'FontName','Helvetica','FontSize',11,'FontWeight','Bold', 'LineWidth', 2);

%perform clustering
eucD = pdist(GroupRF','euclidean');
clustTreeEuc = linkage(eucD,'average');
figure
[~,T3] = dendrogram(clustTreeEuc,2,'Orientation','left','ColorThreshold','default');
clusterIDcommon = T3;

%figure 4f
figure, subplot(2,2,1)
title(strcat('mean RF: ', num2str(mean(GroupRF(clusterIDcommon==1)))))
distributionPlot(percentChange(clusterIDcommon==1,:),'color',[1 1 1],'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
distributionPlot(percentChange(clusterIDcommon==1,:),'xValues',[1.3 2.3 3.3 4.3],'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')

xlabel('tACS frequency')
set(gca, 'XTick',[0.8 2 3.2 4.4])
xticklabels({'0.8' ,'2', '3.2', '4.4'})
hold on
plot(percentChange(clusterIDcommon==1,:)','--','color',[0.6 0.6 0.6],'LineWidth',1)
plot(mean(percentChange(clusterIDcommon==1,:))','g','LineWidth',2)
ylabel('Behavioral Entrainment')
ylim([-5 8])
subplot(2,2,3)
plotSpread(GroupRF(clusterIDcommon==1)')
ylim([0 5])

subplot(2,2,2)
title(strcat('mean RF: ', num2str(mean(GroupRF(clusterIDcommon==2)))))
distributionPlot(percentChange(clusterIDcommon==2,:),'color',[1 1 1],'showMM',2,'addSpread',1,'addBox',1,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
distributionPlot(percentChange(clusterIDcommon==2,:),'xValues',[1.3 2.3 3.3 4.3],'showMM',0,'addSpread',0,'addBox',0,'histOri','right','distWidth',0.2)%,'xyOri','flipped')
xlabel('tACS frequency')
set(gca, 'XTick',[0.8 2 3.2 4.4])
xticklabels({'0.8' ,'2', '3.2', '4.4'})
hold on
plot(percentChange(clusterIDcommon==2,:)','--','color',[0.6 0.6 0.6],'LineWidth',1)
plot(mean(percentChange(clusterIDcommon==2,:))','r','LineWidth',2)
ylabel('Behavioral Entrainment')
ylim([-5 8])
subplot(2,2,4)
plotSpread(GroupRF(clusterIDcommon==2)')
ylim([0 5])