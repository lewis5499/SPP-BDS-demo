function [] = blhpos3d(x,y,z,name,refpos)
%三维散点着色图
clf;

len=length(x);
[bref,lref,href]=xyz2blh(refpos(1),refpos(2),refpos(3));
b=zeros(len,1);
l=zeros(len,1);
h=zeros(len,1);
for i=1:len
    [b(i),l(i),h(i)]=xyz2blh(x(i),y(i),z(i));
end
c=h;%c表示对z轴进行着色

scatter3(b,l,h,50,c,'.');%50表示点的大小，c表示着色情况，'.'表示点的形状
hold on
scatter3(bref,lref,href,300,'red','.');

xlabel('latitude(°)')
ylabel('longitude(°)')
title('3dPos in BLH coordinate system')
grid on
h = colorbar;%右侧颜色栏
set(get(h,'label'),'string','height(m)');%给右侧颜色栏命名

cd ..\imgDir\
saveas(gcf, name, 'png');
cd ..\codeDir\
hold off

end

%xlim([22.4931 23.8784]) %X,Y轴取值范围
%ylim([112.6833 114.5130])

% 
% %导入世界地图
% ax = worldmap('World'); 
% worldmap([-60 60],[0 360])%纬度经度范围显示
% 
% setm(ax, 'Origin', [0 0 0]);
% land = shaperead('landareas', 'UseGeoCoords', true); 
% geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5]);
% % rivers = shaperead('worldrivers', 'UseGeoCoords', true);
% % geoshow(rivers, 'Color', 'blue')
% % cities = shaperead('worldcities', 'UseGeoCoords', true);
% % geoshow(cities, 'Marker', '.', 'Color', 'red')
% %按照经纬度绘制点位
% lat_long = [long lat];
% scatterm(lat_long(:,1),lat_long(:,2),6,'filled','MarkerFaceColor','r')  %经纬度可以是单个点的，也可以是若干个点的
% %plotm(lat,lon,'Marker','.')   %使用scatterm和plotm均可绘制点位



