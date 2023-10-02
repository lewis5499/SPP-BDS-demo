function [] = histogram3(ncats,Y,name)
clf;

X1 = categorical({ncats(1,:),ncats(2,:),ncats(3,:)});
X1 = reordercats(X1,{ncats(1,:),ncats(2,:),ncats(3,:)});
X2 = categorical({ncats(4,:),ncats(5,:),ncats(6,:)});
X2 = reordercats(X2,{ncats(4,:),ncats(5,:),ncats(6,:)});

subplot(1,2,1)
b1=bar(X1,Y(:,1:3));hold on;
b1(1).FaceColor=[0.28 0.57 0.54];
% xtips1 = b1(1).XEndPoints;
% ytips1 = b1(1).YEndPoints;
% labels1 = string(b1(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');
b1(2).FaceColor=[0.73 0.47 0.58];
% xtips2 = b1(2).XEndPoints;
% ytips2 = b1(2).YEndPoints;
% labels2 = string(b1(2).YData);
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');
legend('Single Epoch Least Squares','Kalman Filtering','Orientation','vertical');
xlabel('Categories')
ylabel('Values')
title('Comparison of Kalman Filtering and Single Epoch Least Squares')

subplot(1,2,2)
b2=bar(X2,Y(:,4:6));hold on;
b2(1).FaceColor=[0.28 0.57 0.54];
% xtips1 = b2(1).XEndPoints;
% ytips1 = b2(1).YEndPoints;
% labels1 = string(b2(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');
b2(2).FaceColor=[0.73 0.47 0.58];
% xtips2 = b2(2).XEndPoints;
% ytips2 = b2(2).YEndPoints;
% labels2 = string(b2(2).YData);
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');
legend('Single Epoch Least Squares','Kalman Filtering','Orientation','vertical');
xlabel('Categories')
ylabel('Values')
title('Comparison of Kalman Filtering and Single Epoch Least Squares')

cd ..\imgDir\
saveas(gcf, name, 'png');
cd ..\codeDir\
hold off
end

