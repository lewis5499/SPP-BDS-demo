function [] = RclockdiffDraw(clockmat,name)
 
clf;
set(gcf,'Position',get(0,'ScreenSize'))
% clockmat(:,2)=-1*clockmat(:,2);
%% 连续时间不需处理中断
plot(clockmat(:,1),clockmat(:,2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
scatter(clockmat(:,1),clockmat(:,2),8,[0.28 0.57 0.54],"filled");hold on

%% 非连续时间处理中断
% mark1=[];
% for i=1:length(clockmat)-1
%     if clockmat(i+1,1)-clockmat(i,1)>30
%         mark1=[mark1,i];
%     end
% end
% 
% if isempty(mark1)
%     plot(clockmat(:,1),clockmat(:,2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
% elseif length(mark1)>=2
%     len=length(mark1);
%     plot(clockmat(1:mark1(1),1),clockmat(1:mark1(1),2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
%     for i=1:len-1
%         plot(clockmat(mark1(i)+1:mark1(i+1),1),clockmat(mark1(i)+1:mark1(i+1),2),'Color',[0.28 0.57 0.54],LineWidth=1);
%     end
%     plot(clockmat(mark1(len)+1:end,1),clockmat(mark1(len)+1:end,2),'Color',[0.28 0.57 0.54],LineWidth=1);
% else
%     plot(clockmat(1:mark1(1),1),clockmat(1:mark1(1),2),'Color',[0.28 0.57 0.54],LineWidth=1);hold on
%     plot(clockmat(mark1(1)+1:end,1),clockmat(mark1(1)+1:end,2),'Color',[0.28 0.57 0.54],LineWidth=1);
% end

legend('Receiver clock difference');

% 创建 ylabel
ylabel({'Receiver clock difference (m)'});

% 创建 xlabel
xlabel({'t (sec)'});

% 创建 title
title({name+ ' - t In Bds Sec'});
cd ..\imgDir\
saveas(gcf, name, 'png');
cd ..\codeDir\
hold off
end

