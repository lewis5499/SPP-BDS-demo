function [] = xyzDraw(xyzPosmat,name)

clf;
set(gcf,'Position',get(0,'ScreenSize'))

%% 连续时间不需处理中断
subplot(3,1,1)
plot(xyzPosmat(:,1),xyzPosmat(:,2),'Color',[0.28 0.57 0.54],LineWidth=1,DisplayName='X');hold on
scatter(xyzPosmat(:,1),xyzPosmat(:,2),8,[0.28 0.57 0.54],"filled");
ylabel({'X (m)'});
xlabel({'t (sec)'});
%legend('$X$','Interpreter','latex');
legend('X');
title('Subplot 1: X - t In Bds Sec')

subplot(3,1,2)
plot(xyzPosmat(:,1),xyzPosmat(:,3),'Color',[0.73 0.47 0.58],LineWidth=1,DisplayName='Y');hold on
scatter(xyzPosmat(:,1),xyzPosmat(:,3),8,[0.73 0.47 0.58],"filled");
ylabel({'Y (m)'});
xlabel({'t (sec)'});
legend('Y');
title('Subplot 2: Y - t In Bds Sec')

subplot(3,1,3)
plot(xyzPosmat(:,1),xyzPosmat(:,4),'Color',[0.26 0.45 0.77],LineWidth=1,DisplayName='Z');hold on
scatter(xyzPosmat(:,1),xyzPosmat(:,4),8,[0.26 0.45 0.77],"filled");
ylabel({'Z (m)'});
xlabel({'t (sec)'});
legend('Z');
title('Subplot 3: Z - t In Bds Sec')

cd ..\imgDir\
saveas(gcf, name, 'png');
cd ..\codeDir\
hold off
end


%% 非连续时间处理中断
% mark1=[];
% for i=1:length(xyzPosmat)-1
%     if xyzPosmat(i+1,1)-xyzPosmat(i,1)>30
%         mark1=[mark1,i];
%     end
% end
% 
% %处理中断
% if isempty(mark1)
%     plot(xyzPosmat(:,1),xyzPosmat(:,2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
%     plot(xyzPosmat(:,1),xyzPosmat(:,3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(xyzPosmat(:,1),xyzPosmat(:,4),'Color',[0.26 0.45 0.77],LineWidth=1);
% elseif length(mark1)>=2
%     len=length(mark1);
%     plot(xyzPosmat(1:mark1(1),1),xyzPosmat(1:mark1(1),2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
%     plot(xyzPosmat(1:mark1(1),1),xyzPosmat(1:mark1(1),3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(xyzPosmat(1:mark1(1),1),xyzPosmat(1:mark1(1),4),'Color',[0.26 0.45 0.77],LineWidth=1);
%     for i=1:len-1
%         plot(xyzPosmat(mark1(i)+1:mark1(i+1),1),xyzPosmat(mark1(i)+1:mark1(i+1),2),'Color',[0.28 0.57 0.54],LineWidth=1);
%         plot(xyzPosmat(mark1(i)+1:mark1(i+1),1),xyzPosmat(mark1(i)+1:mark1(i+1),3),'Color',[0.73 0.47 0.58],LineWidth=1);
%         plot(xyzPosmat(mark1(i)+1:mark1(i+1),1),xyzPosmat(mark1(i)+1:mark1(i+1),4),'Color',[0.26 0.45 0.77],LineWidth=1);
%     end
%     plot(xyzPosmat(mark1(len)+1:end,1),xyzPosmat(mark1(len)+1:end,2),'Color',[0.28 0.57 0.54],LineWidth=1);
%     plot(xyzPosmat(mark1(len)+1:end,1),xyzPosmat(mark1(len)+1:end,3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(xyzPosmat(mark1(len)+1:end,1),xyzPosmat(mark1(len)+1:end,4),'Color',[0.26 0.45 0.77],LineWidth=1);
% else
%     plot(xyzPosmat(1:mark1(1),1),xyzPosmat(1:mark1(1),2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
%     plot(xyzPosmat(1:mark1(1),1),xyzPosmat(1:mark1(1),3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(xyzPosmat(1:mark1(1),1),xyzPosmat(1:mark1(1),4),'Color',[0.26 0.45 0.77],LineWidth=1);
%     plot(xyzPosmat(mark1(1)+1:end,1),xyzPosmat(mark1(1)+1:end,2),'Color',[0.28 0.57 0.54],LineWidth=1);
%     plot(xyzPosmat(mark1(1)+1:end,1),xyzPosmat(mark1(1)+1:end,3),'Color',[0.73 0.47 0.58],LineWidth=1);
%     plot(xyzPosmat(mark1(1)+1:end,1),xyzPosmat(mark1(1)+1:end,4),'Color',[0.26 0.45 0.77],LineWidth=1);
% end
% 
% legend('X','Y','Z');
% 
% % 创建 ylabel
% ylabel({'X Y Z (m)'});
% 
% % 创建 xlabel
% xlabel({'t (sec)'});
% 
% % 创建 title
% title({name+' - t In Bds Sec'});



