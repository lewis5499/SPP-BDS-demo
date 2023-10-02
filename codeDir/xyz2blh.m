function [B,L,H] = xyz2blh(X,Y,Z)

%function:
%   将坐标由空间直角坐标系转换为地心坐标系
% input:
%   X、Y、Z-->空间直角坐标系坐标
% output:
%   B、L、H-->地心坐标系坐标

%WGS84
a=6378137;
e2=0.0066943799013;   %第一偏心率的平方

%CGCS2000
%a=6378137;
%e2=0.00669438002290;

%1.计算大地经度L
if (X==0) && (Y>0)
    L=90;
elseif X==0 && Y<0
    L=-90;
else 
    L = atan2(Y,X);
end

%2.利用迭代法计算大地纬度B
%1).赋初值
t_B0=Z/sqrt(X^2+Y^2);           %tanB
t_B1=(a*e2*t_B0/sqrt(1+t_B0^2-e2*t_B0^2)+Z)/sqrt(X^2+Y^2);

%2).进行迭代
while abs(t_B1-t_B0)>1e-10
    t_B0=t_B1;
    t_B1=(a*e2*t_B0/sqrt(1+t_B0^2-e2*t_B0^2)+Z)/sqrt(X^2+Y^2);
end
B=atan2(t_B1,1);

%3).计算大地高
N=a/sqrt(1-e2*(sin(B))^2);
H=sqrt(X^2+Y^2)/cos(B)-N;

