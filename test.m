%close all;
clc;

x = [0;0.25;0.5;0.75;1];
y = [0.1;0.5;0.75;0.6;-0.1];
a = min(x);
b = max(x);
np = numel(x);

nCurves = 3;
curveOrder = 1;
points = [y,x];
[fhandle,exitflag,message] = curvefit(points,nCurves,curveOrder,'plot',true,'floating',true,'beta',40);

x2 = linspace(a,b,max(np*10,1000))';
fhandle(x2,'plot',true,'df',[1,2,3]);

return
clear all;
close all;
clc;
np = 15;
a = pi/2;
b = 2*pi;
x = linspace(a,b,np)';
y = sin(x);

nCurves = 5;
curveOrder = 4;
points = [y,x];
pointlb = y*1.05;
pointub = y-y*0.2;
c0 = [a,0; b,0];
c1 = [a,0; b,0];
[fhandle,exitflag,message] = curvefit(points,nCurves,curveOrder,'plot',true,'floating',true,'beta',100);

x2 = linspace(a,b,max(np*10,1000))';
fhandle(x2,'plot',true,'df',[1,2,3]);