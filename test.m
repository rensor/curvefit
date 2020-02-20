close all;
clc;

x = [0;0.25;0.5;0.75;1];
y = [0.1;0.5;0.75;0.6;-0.1];

nCurves = 5;
curveOrder = 1;
points = [y,x];
[fhandle,exitflag,message] = curvefit(points,nCurves,curveOrder,'plot',true);



return
%clear all;
close all;
clc;
np = 20;
x = linspace(0,2*pi,np)';
y = sin(x);

nCurves = 10;
curveOrder = 3;
points = [y,x];
pointlb = y*1.05;
pointub = y-y*0.2;
c0 = [0,0; 2*pi,0];
c1 = [0,0; 2*pi,0];
[fhandle,exitflag,message] = curvefit(points,nCurves,curveOrder,'plot',true);

x2 = linspace(0,2*pi,np*10)';
fhandle(x2,'plot',true,'df',[1,2,3]);