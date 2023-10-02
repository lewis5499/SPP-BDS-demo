function [] = m0Draw(m0,name)

clf;
set(gcf,'Position',get(0,'ScreenSize'))

%% 连续时间观测无需处理中断
plot(m0(:,1),m0(:,2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
plot(m0(:,1),m0(:,3),'Color',[0.73 0.47 0.58],LineWidth=1);

scatter(m0(:,1),m0(:,2),8,[0.28 0.57 0.54],"filled");hold on
scatter(m0(:,1),m0(:,3),8,[0.73 0.47 0.58],"filled")

legend('singleLS','staticKF');
ylabel({'m0 (m)'});
xlabel({'t (sec)'});
title({name+' - t In Bds Sec'});

cd ..\imgDir\
saveas(gcf, name, 'png');
cd ..\codeDir\
hold off
end

%% 不连续时间观测处理中断
% mark1=[];
% for i=1:length(neuPosmat)-1
%     if neuPosmat(i+1,1)-neuPosmat(i,1)>30
%         mark1=[mark1,i];
%     end
% end
%
% if isempty(mark1)
%     plot(neuPosmat(:,1),neuPosmat(:,2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
%     plot(neuPosmat(:,1),neuPosmat(:,3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(neuPosmat(:,1),neuPosmat(:,4),'Color',[0.26 0.45 0.77],LineWidth=1);
% elseif length(mark1)>=2
%     len=length(mark1);
%     plot(neuPosmat(1:mark1(1),1),neuPosmat(1:mark1(1),2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
%     plot(neuPosmat(1:mark1(1),1),neuPosmat(1:mark1(1),3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(neuPosmat(1:mark1(1),1),neuPosmat(1:mark1(1),4),'Color',[0.26 0.45 0.77],LineWidth=1);
%     for i=1:len-1
%         plot(neuPosmat(mark1(i)+1:mark1(i+1),1),neuPosmat(mark1(i)+1:mark1(i+1),2),'Color',[0.28 0.57 0.54],LineWidth=1);
%         plot(neuPosmat(mark1(i)+1:mark1(i+1),1),neuPosmat(mark1(i)+1:mark1(i+1),3),'Color',[0.73 0.47 0.58],LineWidth=1);
%         plot(neuPosmat(mark1(i)+1:mark1(i+1),1),neuPosmat(mark1(i)+1:mark1(i+1),4),'Color',[0.26 0.45 0.77],LineWidth=1);
%     end
%     plot(neuPosmat(mark1(len)+1:end,1),neuPosmat(mark1(len)+1:end,2),'Color',[0.28 0.57 0.54],LineWidth=1);
%     plot(neuPosmat(mark1(len)+1:end,1),neuPosmat(mark1(len)+1:end,3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(neuPosmat(mark1(len)+1:end,1),neuPosmat(mark1(len)+1:end,4),'Color',[0.26 0.45 0.77],LineWidth=1);
% else
%     plot(neuPosmat(1:mark1(1),1),neuPosmat(1:mark1(1),2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
%     plot(neuPosmat(1:mark1(1),1),neuPosmat(1:mark1(1),3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(neuPosmat(1:mark1(1),1),neuPosmat(1:mark1(1),4),'Color',[0.26 0.45 0.77],LineWidth=1);
%     plot(neuPosmat(mark1(1)+1:end,1),neuPosmat(mark1(1)+1:end,2),'Color',[0.28 0.57 0.54],LineWidth=1);
%     plot(neuPosmat(mark1(1)+1:end,1),neuPosmat(mark1(1)+1:end,3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(neuPosmat(mark1(1)+1:end,1),neuPosmat(mark1(1)+1:end,4),'Color',[0.26 0.45 0.77],LineWidth=1);
% end