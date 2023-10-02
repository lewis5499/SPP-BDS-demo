function pos_NEU = xyz2neu(p0_XYZ,p1_XYZ)

%function:
%   将p1点坐标由空间直角坐标系转换为站心坐标系(站心坐标系原点为p)
% input:
%   p0_XYZ-->测站点坐标
%   p1_XYZ-->空间点坐标
% output:
%   N、E、U-->地心坐标系坐标

% 1.构造旋转矩阵R
x0=p0_XYZ(1); y0=p0_XYZ(2); z0=p0_XYZ(3);
x1=p1_XYZ(1); y1=p1_XYZ(2); z1=p1_XYZ(3);

%调用函数得到纬度B、经度L
[B,L,~]=xyz2blh(x0,y0,z0);        

R=[-sin(B)*cos(L) -sin(B)*sin(L)  cos(B)
      -sin(L)         cos(L)        0
   cos(B)*cos(L)   cos(B)*sin(L)  sin(B)];

% 2.实现转换
pos_NEU=(R*[x1-x0;y1-y0;z1-z0]).';

end
