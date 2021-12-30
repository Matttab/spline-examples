% Test of some routines from linuxCNC source code to interpolate splines
% Check of spline construction integrity
% Spline from P0 to P3 with interior points P1,P2 as waypoint/node/knot
%
% SPLINE_FEED with 6 input arguments require as input
% the control points of a B-Spline. P0 and P3 are as well
% use as first and last control points, however, only
% [P1,P2,P3] are passed to SPLINE_FEED as arguments.

% P1 and P2 define the shape of the spline
% There are four shape parameters P1x, P1y, P2x, P2y
% Likely the x-grid is always uniform from 
% P0123 = P0+[0 1/3 2/3 1]*(P3-P0)];

% x0,y0 Start point, stored in class property _pos_x and _pos_y
%
% Spline defition in Bezier format. Arguments of SPLINE_FEED 
% y1 first coeff of cubic B-Spline
% y2 second coeff of cubic B-Spline
% y3 third coeff of cubic B-Spline

% x0,y0 first control point on outer polygon = Start point P0
% x1,y1 second control point on outer polygon = P1
% x2,y2 third control point on outer polygon  = P2
% x3,y3 last control point on outer polygon = End point P3

clear all; close all;
mySpline = linuxCNC;

mySpline.metric=true;
x0=0;
y0=0;
mySpline._pos_x=x0;
mySpline._pos_y=y0;

% P123 =  [I, P, endpt_x; J, Q, endpt_y]
P123  = [1.0, 2.0, 3.0,;
         1.0, 1.7, 1.8]; 
P0123=[[x0; y0] P123];

% scale it to x=[0..1]
%P123  = [11/3, 2/3, 1;
%         1.0, 1.7, 1.8]; 
%P0123=[[x0; y0] P123];


% do a spline interpolation separately for x and y
t=linspace(0,1,length(P0123));
tt=linspace(0,1,length(P0123)*10+1);
ppx=splinefit(t,P0123(1,:), t ,"order", 3); 
ppy=splinefit(t,P0123(2,:), t ,"order", 3); 
xxp=ppval(ppx,  tt);
yyp=ppval(ppy,  tt);

% use x as parameter/breaks and fit a spline to y data
pp=splinefit(P0123(1,:), P0123(2,:),t ,"order", 3); 
xx=tt*(P0123(1,end)-P0123(1,1))+P0123(1,1);
yy=ppval(pp,xx);
% from octave
% pp =  [ 1.69  0.38 0.01 0.0; % first segment is almost linear
%        -2.60  2.08 0.91 1.37;
%         0.05 -0.52 1.44 0.58]

% matlab comparison
% qq= spline(P0123(1,:), P0123(2,:));
% also with breaks = [0 1 2 3 4];
% will give other coefficients than octave splinefit
qq.coefs =  [-0.05  0.00 1.05 0.0; % first segment is almost linear
        -0.05 -0.15 0.90 1.0;
        -0.05 -0.30 0.45 1.7];

% Conversion to B-Spline using Matlab
% bpp = fn2fm(qq, 'B-', 0.1)
% will change the breaks to knots with 4-fold multiplicity at start and end pt
% bpp.knots = [0 0 0 0 3 3 3 3]
bpp.coefs = [0 1.05 2.1 1.8] ; %the corners of the control polygon y-coordinates
[xb,yb]=mySpline.SPLINE_FEED3(1, bpp.coefs(2), 2, bpp.coefs(3), 3, bpp.coefs(4) ); 

% For a set of G5 codes
% Several splines in a row        
[x4, y4]= mySpline.construct_spline(P0123);

% We want only one spline segment between start and endpt
[x, y]= mySpline.construct_spline([P0123(:,1) P0123(:,4)]);
  
         
ctrlPt(:,1)=P0123(:,1);       
ctrlPt(:,2)=[1; 2];   
ctrlPt(:,3)=[1.05; 2.1];   
ctrlPt(:,4)=P0123(:,4);
         
[xq,yq]=mySpline.SPLINE_FEED2(ctrlPt(1,2), ctrlPt(2,2),   ctrlPt(1,4),  ctrlPt(2,4)  );
[xc,yc]=mySpline.SPLINE_FEED3(ctrlPt(1,2), ctrlPt(2,2),   ctrlPt(1,3),  ctrlPt(2,3),  ctrlPt(1,4),ctrlPt(2,4) ); 

figure; hold on;
plot(P0123(1,:), P0123(2,:),'-o');

%plot(xq,yq,'--r');
plot(x,y,':k');
plot(xb,yb,'-c');
%plot(xc,yc,'-g');
plot(xx,yy,'-m');
%plot(xxp,yyp,'-r');
plot(P0123(1,:),bpp.coefs ,'-dc');

for k=1:3
 plot(x4(k,:),y4(k,:),':k');
end
legend({'nodes',  ...
'spline constructor (output of 1 splines from start to end pt)', ...
'SPLINE FEED3(x of nodes, coefs from MATLAB converson to B-Spline)', ...
'cubic octave fit explicit y=y(x)', ...
'Coefs=control polygon from matlab conversion to B-Spline', 'spline constructor (output of 3 splines)',}, 'location','southeast'  );

% conversion works only on
% parameterspace t=[0..1]
% transformation x=x(t)
% to be implemented
linuxCNC.kanonical2bernstein(qq.coefs)
linuxCNC.bernstein2kanonical(ctrlPt)





