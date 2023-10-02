function [] = histogram6(ncats,Y,name)
clf;

X = categorical({ncats(1,:),ncats(2,:),ncats(3,:),ncats(4,:),ncats(5,:),ncats(6,:)});
X = reordercats(X,{ncats(1,:),ncats(2,:),ncats(3,:),ncats(4,:),ncats(5,:),ncats(6,:)});

b=bar(X,Y);
b(1).FaceColor=[0.28 0.57 0.54];
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
b(2).FaceColor=[0.73 0.47 0.58];
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
legend('Single Epoch Least Squares','Kalman Filtering','Orientation','vertical','Location','northwest');

xlabel('Categories')
ylabel('Values')
title('Comparison of Kalman Filtering and Single Epoch Least Squares')
cd ..\imgDir\
saveas(gcf, name, 'png');
cd ..\codeDir\
hold off
end

