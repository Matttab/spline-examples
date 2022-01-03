% Test of FEED routine from linuxCNC source to interpolate splines
% ----------------------------------------------------------------------
% SPLINE_FEED(P1x, P1y, P2x, P2y, P3x, P3y) with inputs
% P1 and P2: interior control points of cubic spline
% P3: end point of spline feed.
% P0=[x0,y0]: start point stored in class properties [_pos_x,_pos_y]
% P1-P0 defines slope at start point as (y1-y0)/(x1-x0) at start
% P3-P2 defines slope at end point as (y3-y2)/(x3-x2) at end
% ----------------------------------------------------------------------
% SPLINE_FEED(P1x, P1y, P2x, P2y) with inputs
% P1: interior control point of quadratic spline
% P2: end point of spline feed.
% P0=[x0,y0]: start point, stored in class property [_pos_x,_pos_y]
% P1-P0 defines slope at start point as (y1-y0)/(x1-x0) at start
% P2-P1 defines slope at end point  as (y2-y1)/(x2-x1) at end
% ----------------------------------------------------------------------
% P1 and P2 define the shape of the spline
% There are four shape parameters P1x, P1y, P2x, P2y for cubic spline
% The x-grid can but does not have to be uniform, e.g. 
% P0123 = P0+[0 1/3 2/3 1]*(P3-P0)] can be used.
% ----------------------------------------------------------------------
% Spline defition is in Bezier format. Arguments of SPLINE_FEED 
% y1 first coeff of cubic B-Spline
% y2 second coeff of cubic B-Spline
% y3 third coeff of cubic B-Spline
% ----------------------------------------------------------------------
% P0=[x0,y0] first control point of cubic B-Spline
% P1=[x1,y1] first control point of cubic B-Spline
% P2=[x2,y2] first control point of cubic B-Spline
% P3=[x3,y3] first control point of cubic B-Spline
% ----------------------------------------------------------------------
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

% Create double multiplicity for start and end points
% to add constraints for y'(x)=dy/dx at start and end
% 
con.xc = P0123(1,[1,1,end,end]);% constraint locations
con.yc = [1.05 P0123(2,[1,end]) -0.3]; % constraint values for y,y',y''
con.cc = [0 1 1 0;1 0 0 1]; %constrain [y(start) 'y(start) y(end) y'(end)]

tt=linspace(0,1,length(P0123)*10+1); % uniform interpolation grid

% use x as parameter/breaks and fit a spline to y data
% regrid to uniform grid
pp=splinefit(con.xc, con.yc,con.xc ,"order", 3,'constraints',con); 
xx=tt*(P0123(1,end)-P0123(1,1))+P0123(1,1);
yy=ppval(pp,xx);
yy(1);
yy(end);
dd=diff(yy)./diff(xx);
%dd(1) = 1.0497 converges to 1.05 for dense [0 1] grid
%dd(end)=-0.26653 converges to -0.3for dense [0 1] grid
% from splinefit including proper boundary constraints
% pp =  [-0.05 0 1.05 0];

% matlab comparison
% qq= spline(P0123(1,:), P0123(2,:));
% with breaks at t= [0 1 2 3];
% will give other coefficients than octave splinefit
qq.coefs =  [-0.05  0.00 1.05 0.0; % first segment is almost linear
        -     0.05 -0.15 0.90 1.0;
             -0.05 -0.30 0.45 1.7];

% Conversion to B-Spline using Matlab
% bpp = fn2fm(qq, 'B-', 0.1)
% will change the breaks to knots with 4-fold multiplicity at start and end pt
% bpp.knots = [0 0 0 0 3 3 3 3]
bpp.coefs   =[0 1.05 2.1 1.8] ; % the corners of the control polygon y-coordinates
% for      x=[0 1.00 2.0 3.0]   % the diff gives the slope at start and end points
[xb,yb]=mySpline.SPLINE_FEED3(1, bpp.coefs(2), 2, bpp.coefs(3), 3, bpp.coefs(4) ); 
db=diff(yb)./diff(xb);
%db(1) = 1.05 = (y1-y0)/(x1-x0) at start
%db(end)=-0.3 = (y3-y2)/(x3-x2) at end

% For a set of G5 codes
% Several splines in a row  
% This constructs a convex superposition
% but does not respect the boundary conditions!     
[x4, y4]= mySpline.construct_spline(P0123);

% We acwant only one spline segment between start and endpt
[x, y]= mySpline.construct_spline([P0123(:,1) P0123(:,4)]);
  
         
ctrl3Pt(:,1)=P0123(:,1);       
ctrl3Pt(:,2)=[1; 1.05];   
ctrl3Pt(:,3)=[2; 2.1];   
ctrl3Pt(:,4)=P0123(:,4);
%diff(ctrl3Pt(2,:))./diff(ctrl3Pt(1,:))
ctrl2Pt(:,1)=P0123(:,1);       
ctrl2Pt(:,2)=[2; 2.1];     
ctrl2Pt(:,3)=P0123(:,4);
%diff(ctrl2Pt(2,:))./diff(ctrl2Pt(1,:))      
[xq,yq]=mySpline.SPLINE_FEED2(ctrl2Pt(1,2), ctrl2Pt(2,2),   ctrl2Pt(1,3),  ctrl2Pt(2,3)  );
[xc,yc]=mySpline.SPLINE_FEED3(ctrl3Pt(1,2), ctrl3Pt(2,2),   ctrl3Pt(1,3),  ctrl3Pt(2,3),  ctrl3Pt(1,4), ctrl3Pt(2,4) ); 

figure; hold on;
plot(P0123(1,:), P0123(2,:),'+');

plot(xq,yq,'--r');
%plot(x,y,':k');
%plot(xb,yb,'-b');

plot(xx,yy,'-m');
plot(xc,yc,'--b');
%plot(xxp,yyp,'-r');
plot(P0123(1,:),bpp.coefs ,'-db');

for k=1:3
 plot(x4(k,:),y4(k,:),'--k');
end
legend({'nodes',  ...
"SPLINE FEED2(y\'(0)=1.05,y\'(1)=-0.3  )", ...
'cubic octave fit explicit y=y(x)', ...
"SPLINE FEED3(y\'(0)=1.05,y\'(1)=-0.3  )", ...
'Coefs=control polygon for B-Spline', 'spline constructor (convex sum of 3 splines)',}, 'location','southeast'  );

% conversion works for
% parameterspace t=[0..1]
linuxCNC.kanonical2bernstein(pp.coefs)
linuxCNC.bernstein2kanonical([3.1 -6.3 3.15 0])





