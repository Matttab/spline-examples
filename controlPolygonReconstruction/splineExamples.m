% Examples to work with splines in MATLAB
% Kind of sandbox to experiment with different spline types
% 
% Fit a piecewise cubic spline (of type pp) with breaks
% pp = spline (x, y);  using x as parameter for y=y(x)
% yy = ppval(pp, xx);  interpolate y to other values xx
%  output: pp.coefs, pp.breaks
%  make a pp spline from known coefs: ppmak(breaks,coefs) 
% ----------------------------------------------------------------------
% Use of MATLAB 'curve fitting' toolbox

% pp = csapi(parameter,[x;y]) make a pp spline
%   use csape() for periodic bc=boundary conditions
%  
% sp = fn2fm(pp,'B-',0.1); %conversion to B-Spline
% pp = fn2fm(sp,'pp',0.1); %conversion to pp Spline
%
% sp = spmak(knots, coefs) make a B-Spline from known coefs
%   output: sp.coefs, sp.knots (knots instead of breaks)
%
% sp = spmak(parameter,[x;y])
%   output: sp.coefs, sp.knots
%
% Challenge: Reconstruct the control node polygon of a given pp spline.
% Use this polygon to construct a B-Spline similar to the inital pp spline.
%
% To be improved: Padding of nodes to define boundary condition
%  Any machine based approximation will require an automated way
%  to define boundary conditions for an arbitrary shape correctly.
% Manual padding additional nodes?
%  Challenge: For which shape to add nodes where exactly ?
%
% The code of the published source
% https://www.mathworks.com/matlabcentral/answers/86615-how-to-plot-an-ellipse#answer_96122
% has been used as starting base.
%%
% Example 1: some odd peak shape to be fit by spline
% Parts of code are reprint from The MathWorks Inc. help center
% ----------------------------------------------------------------------
clear all; close all;
% Create waypoints to interpolate with a cubic Spline.
% how to define the boundary slope: [y' y y']
% derivates with respect to parameter y=y(t), y'=dy/dt
wptsx = [-5 -4:4 5];
wptsy = [0 0 0 1.12 2.36 2.36 1.46 .49 .06 0.06 0.06];

% for this curve we can use x as parameter
% this is not possible for multple occurance of same x value
% y must be uniquely defined for each parameter value
tt=linspace(min(wptsx),max(wptsx),31);
%define slopes at bounary as 0

pp1= spline(wptsx,[0 wptsy 0]); 
pp2= spline(wptsx,wptsy);

% spline command is same as csapi
% pp3= csapi(wptsx,wptsy);
% pp4 = csape( wptsx, wptsy, 'periodic' );
% s1 = spline(wptsx,[0 wptsy 0],tt);

% interpolate to higher resolution
yy1 = ppval(pp1, tt);
yy2 = ppval(pp2, tt);

% plot
figure; hold on;
    plot(tt,yy1,'-sk');
    plot(tt,yy2,'-dg');
    plot(wptsx,wptsy,'-ob');

% check the slopes at boundaries
slope1=diff(yy1)./diff(tt);
slope2=diff(yy2)./diff(tt);
slope1([1,end])
slope2([1,end])
abs(slope1-slope2)./max(abs(slope2),0)

figure; hold on;
    plot(tt(1:end-1)+0.5*(tt(2)-tt(1)),slope1,'-sk')
    plot(tt(1:end-1)+0.5*(tt(2)-tt(1)),slope2,'-dg')

%%    
% Example 2: Approximate a circle with cubic Splines
% Parts of code are reprint from The MathWorks Inc. help center
% ----------------------------------------------------------------------
phi = pi*[0:.5:2]; 
pts = [0  1  0 -1  0  1  0;  %x
       1  0  1  0 -1  0  1]; %y
   % two more column to fix the bc x'=0 and y'=1
   % only 5 = length(phi) points will be used to
   % cover [0..2pi]
   
% the proper control nodes for a B-Spline are more tricky,
% in case you want to reconstruct from number=7 pts use
pts7=1.5*[0.6667    0.6667    0.0000   -1.0000   -0.0000    0.6667    0.6667;
          0.0000    0.3491    0.9933    0.0000   -0.9933   -0.3491   -0.0000];
      
%Simple approach fo B-Slpine with number=8 pts:
% We use the quadrant (+/-1 0) pts without the
%          appended bc of [x',y'] = [0 1]
% We double the array to fit a double loop of 8 pts
% Used is actually only the part from angle=[0..2pi]
% The parameters needs to cover the full double loop
% thus phi = [-1..3]pi, at each bc 1 pi lead to fix bc
bPhi=pi*linspace(-1,3,12);
bB=spmak(bPhi,1.5*[pts(:,2:5) pts(:,2:5)]); 
% the factor 1.5 is a riddle to be figured out

pp = spline(phi,pts); %separate fit for each column in pts
diff(pts)/(0.5*pi)

phi3 = phi
pp3 = csapi(phi3,pts(:,2:6));  %this will fail due to undefined bc
% use csape for periodic bc
% this works nice and does not need the bc to be fixed
% only pts need are the quadrant pts covering [0..2]*pi
pts3p = [1  0 -1  0  1;  %x
         0  1  0 -1  0]; %y
pp3p = csape( phi3, pts3p, 'periodic' ); 

% circle
phii=linspace(min(phi),max(phi),101);
cc=[cos(phii); sin(phii)];

% plot
yy = ppval(pp, phii);
yy3 = ppval(pp3, phii);
yy3p = ppval(pp3p, phii);
yB= fnval(bB, phii);
figure; hold on;
    plot(cc(1,:),    cc(2,:), ':k');
    plot(yy(1,:),    yy(2,:),'-b')
    plot(yy3(1,:),   yy3(2,:),':m')
    plot(yy3p(1,:),  yy3p(2,:),'--r')
    plot(pts(1,2:5), pts(2,2:5),'or')
    plot(yB(1,:),    yB(2,:),'.g')
    axis equal
    
%%
% Example 3.1: Approx a circle with cubic B-Splines
% Parts of code are reprint from The MathWorks Inc. help center
% ----------------------------------------------------------------------
points = [0 -1 0 1;              %x
          1 0 -1 0]*1.6851;      %y

% Simple but not best approach to keep number of knots low
% We do pad the knot vector with a second 360 degree period
% Parameter range is arbitray, here taken -4:8 to map range of 4*pi.
% Actually good will be the middle 0..4 maping to 0.pi
% There is a lead-in of pi, from -4..0
% and a lead-out of pi,     from 4..8
% Nicer would be a periodic definition with no need for data padding

sp = spmak(-4:8,[points points]);
sp.coefs
sp.knots
%sp = spmak([0:0.1:1]*2*pi,[points points]);
% provides a planar, quartic, spline curve whose middle part is a pretty good approximation to a circle, as the plot on the next page shows. It is generated by a subsequent
figure;
plot(points(1,:),points(2,:),'o'), hold on 
fnplt(sp,[0,4]), axis equal square, hold off

% get min and max curvature
% curvature kappa = abs(x'y''-y'x'')/(x'^2+y'^2)^(3/2)
t = linspace(0,4,21);zt = zeros(size(t));
dsp = fnder(sp); dspt = fnval(dsp,t); ddspt = fnval(fnder(dsp),t);
kappa = abs(dspt(1,:).*ddspt(2,:)-dspt(2,:).*ddspt(1,:))./(sum(dspt.^2)).^(3/2);
[min(kappa),max(kappa)] 
% kappa = [1.6747    1.8611]
% 1/norm(fnval(sp,0))= 1.7864 

%%

% ----------------------------------------------------------------------
% Example 3.2: Conversion of pp to B-spline for case of circle
% Create circle with a cubic B-Spline. Convert to B-spline and back.
% Parts of code are reprint from The MathWorks Inc. help center
% ----------------------------------------------------------------------
% g = fn2fm(f,form)
% describes the same function as is described by f, but in the form specified
% by  'B-', 'pp', 'BB', 'rB', 'rp',
% for the B-form, the ppform, the BBform, and the two rational spline forms.

% Simple but not best approach to keep number of knots low
% We do pad the knot vector with a second 360 degree period

pts = [0  1  0 -1  0  1  0 ;  %x
       1  0  1  0 -1  0  1 ]; %y
   % two more column to fix bc
   
pts2 = [0  1  0 -1  0  1 ;  %x
       1  0  1  0 -1  0]; %y
   % two more column to fix bc

pts3 = [0  1  0 -1  0  1  0 -1;  %x
       1  0  1  0 -1  0  1  0]; %y
   % two more column to fix bc
   
phi = linspace(0,2,length(pts)-2)*pi; 
pp = spline(phi,pts);

% conversion of pp to B-spline
bpp = fn2fm(pp,'B-',0.1)
ppbpp = fn2fm(bpp,'pp',0.1)

phii = linspace(0,2,101)*pi; 
yy = ppval(pp, phii);
byy= fnval(bpp, phii);

% circle
phii=linspace(min(phi),max(phi),101);
cc=[cos(phii); sin(phii)];

%plot
figure; hold on;
    plot(yy(1,:),yy(2,:),'-k')
    plot(-yy(1,:),yy(2,:),'-g')
    plot(byy(1,:),byy(2,:),'-b')
    plot(pts(1,2:5),pts(2,2:5),'or')
    plot(bpp.coefs(1,:),bpp.coefs(2,:),'-xm')
    plot(cc(1,:),    cc(2,:), ':k');
    axis equal

% breaks vs. knots
% and zero padding for B-spline
% first 3 and last 3 elements
% order 4
% pieces 4 for cubic spline
% number 4+3 for B-spline 
pp.breaks
bpp.knots(4:end-3)
ppbpp.breaks

pp.coefs
bpp.coefs
% some error of 1E-15 in coefs after back conversion
(ppbpp.coefs-pp.coefs)/1E-15

% the control points of the B-spline 
% are asymmetrically placed - SURPRISE !
% control points are asymmetric to get a symmetric shape
% bpp.coefs =
% 1.0000    1.0000    0.0000   -1.5000   -0.0000    1.0000    1.0000
%      0    0.5236    1.4899    0.0000   -1.4899   -0.5236   -0.0000

% The circle diameter to control point diameter
% have a ratio of 1.9

% bpp.coefs/1.5 =

% pts2=[    0.6667    0.6667    0.0000   -1.0000   -0.0000    0.6667    0.6667;
%          0    0.3491    0.9933    0.0000   -0.9933   -0.3491   -0.0000];

%%
% Example 4:  Interpolate with B-Spline vs. cubic pp spline
% ----------------------------------------------------------------------
wptsx = [0 1 2.1 8 4 3];
wptsy = [0 1 1.3 .8 .3 .3];
wpts = [wptsx; wptsy];
L = length(wpts) - 1;
%t = 0:0.05:1;

% spline cannot be directly used on (x,y)
% can do it with separate spline for x and y
tt=linspace(0,1,101);
t=linspace(0,1,length(wptsx));

% Show spline interpolation of x and y (done separately)
% figure; hold on;
% plot(t, wptsx,'-k')
% plot(t, wptsy,'-b')

% check boundary slopes
diff(wptsx)./diff(t);
diff(wptsy)./diff(t);
% define slopes at bounary same as at
% point before boundary

% cubic spline
pp1x= csape(t, [5 wptsx -5]);
pp1y= csape(t, [5 wptsy 0]); 
xx1 = ppval(pp1x, tt);
yy1 = ppval(pp1y, tt);
% B-Spline
pp3x = spapi(3,t, wptsx);
pp3y = spapi(3,t, wptsy);
xx3 = fnval(pp3x, tt);
yy3 = fnval(pp3y, tt);

figure; hold on
    plot(xx3, yy3, '-b');
    plot(xx1, yy1, '-k');
    plot(wptsx, wptsy, 'dr');

%%
% ----------------------------------------------------------------------
% Example 5.1:  How to find the outer envelope polygon ?
% There should be a geometric solution to find the control polygon.
% We will work back from a solution extracted from mathworks help
% center. An explanation of this approach or any interpretation is
% however not given in the help.
% Curious to figure it out yourself? Then do not read further.
% ----------------------------------------------------------------------
% Create waypoints to interpolate with a B-Spline.
wptsx = [0 1 2.1 8 6 3];
wptsy = [0 1 1.3 .8 .3 .3];
% Put in vector to solve later a linear equation system LES.
wpts = [wptsx; wptsy];
L = length(wpts) - 1;

% Form matrices used to compute interior points of control polygon
% Interior points are just the normal pts excluding the start and end pt
r = zeros(L+1, size(wpts,1));
A = eye(L+1);
for i= 2:L
    A(i,(i-1):(i+1)) = [1 4 1];  %1+4+1 = 6
    r(i,:) = 6*wpts(:,i)';       %tranposed vector wpts
end
% Somewhat arbitray, the definition of additonal bc
% Override end points and choose r0 and rL.


r(1,:) = (wpts(:,1) + (wpts(:,2) - wpts(:,1))/2)';
r(end,:) = (wpts(:,end-1) + (wpts(:,end) - wpts(:,end-1))/2)';

% Here come the LES = LINEAR EQUATION SYSTEM
% This system can be solved geometrically by drawing lines and midpoints.
% ----------------------------------------------------------------------
% Right division: solve linear equations system LES using e.g. Gauss
% method.
%   A*dInterior = r;
    dInterior = (A\r)';
    wptsInterior=A*dInterior'/6
% Construct a complete control polygon adding
% the intial and final knot. They will serve as well as corners of the polygon:
    cpts = [wpts(:,1) dInterior wpts(:,end)];
% ----------------------------------------------------------------------
% the EQs for x=dInterior read
%   x2 + 4 x3 + x4 = r2y = 6*wptsx3
%   y2 + 4 y3 + y4 = r2y = 6*wptsy3

% known: the midpoints of x2 and x3 should be close to the spline curve
% we try to rewrite the LES to get some interpretation:
%   x2 + 2 x3 + x4 = r2y = 6*wptsx3 - 2 x3 
%   y2 + 2 y3 + y4 = r2y = 6*wptsy3 - 2 y3 

% to be figured out, is the meaning of the relation of the 
% interior polygon points x2,x3,....x(N-1) to each other,
% and with respect to the inital knots.

% A trial to figure out the meaning of the LES
% ----------------------------------------------------------------------
% divide by 2 and group to midpt values
%   0.5*(x2 + x1) + 0.5*(x2+x3) = 2*wptsx3 + (wptsx3 - x3)
%   0.5*(y2 + y1) + 0.5*(y2+y3) = 2*wptsy3 + (wptsy3 - y3)

% divide by 2 and group again to midpt values
%   0.5*[0.5*(x2 + x1) + 0.5*(x2+x3)] - 0.5*(wptsx3 - x3) = wptsx2  
%   0.5*[0.5*(y2 + y1) + 0.5*(y2+y3)] - 0.5*(wptsy3 - y3) = wptsy2  
% The midpt of two subsequent polygon midpoints, are some distance aways from the
% original knots. This distance is half of (wptsx3 - x3),
% the distance of the polygon corner from its respecting knot.
% The 3 pts lie on a line connecting, for each knot, the
% midpt of two neighboring polygon midpoints and polygon corner
% through the knot, intersecting it at ratio 1:2.
% The polygon corner is 3 times the distance from the
% midpt of two neighboring polygon midpoints to the waypts
% away.
% At least this is true for all interior points (unaffected by bc)

% some pts for 1st order spline we can read from our initial constraints.
q1(1,:) =[ 0 1 2.1 8 6 3];
q1(2,:) = [0 1 1.3 .8 .3 .3];

% the other pts for dInterior we should insert such that the existing
% pts q1 close to the midpoints on the polygon edges.

% Helpful is bsplinepolytraj to compute a polynomial with the new control points
% Unfortunately, this function is in the inaccesible "robotic systems"
% toolbox and not as expected, in the "curve fitting" aka spline/nurbs toolbox

% CANNOT be calculated without "robotic system" toolbox
%q = bsplinepolytraj(cpts, [0 1], t); 
% ----------------------------------------------------------------------
% we will use a cubic spline instead to plot something
% that follow the inital waypts = knots
t=linspace(-1/6,7/6,8);
tt=linspace(0,1,101);
% diff([0 wptsx 0])./diff(t);
% diff([0 wptsy 0])./diff(t);
pp1x= spline(t,[-1 wptsx 1]); %probably different bc than above
pp1y= spline(t,[-1 wptsy 0]); %probably different bc than above

xx1 = ppval(pp1x, tt);
yy1 = ppval(pp1y, tt);

% Plot the results. Show the original waypoints, the computed polygon, and the interpolated B-spline.
figure;
hold all
plot(wptsx, wptsy, 'ob');
plot(xx1, yy1, '-k');
plot(cpts(1,:), cpts(2,:), '-sb');
plot(0.5*(cpts(1,1:end-1)+cpts(1,2:end)), 0.5*(cpts(2,1:end-1)+cpts(2,2:end)), 'xb');
plot(q1(1,:), q1(2,:),'-r');
legend('Original waypoints', 'Cubic spline', 'Computed control polygon');
% ----------------------------------------------------------------------
% from the online life editor we can get
% for cpts above, tpts=[0 1], tvec = 0:1/5:1
% and [q2, qd, qdd, pp] = bsplinepolytraj(cpts,tpts,tvec);

q2(1,:) = [0 1 2.1 8     6  3];
q2(2,:) = [0 1 1.3 0.8 0.3 0.3];

q3=[];
for k=1:6
    q3(:,1+(k-1)*2) = q2(:,k);
end
% manually we add the intermediate pts for tvec = 0:1/10:1
q3(1,2:2:10) = [0.6809 0.9832 5.2988 7.7338 4.6532];
q3(2,2:2:10) = [0.6080 1.2474 1.1148 0.4535 0.2731];
% ----------------------------------------------------------------------
%[1] Farin, Section 9.1
% Copyright 2018 The MathWorks Inc.
% ComputeBSplineTrajectoryFor2DPlanarMotionExample.mlx

% Copyright 2020 The MathWorks Inc.
% openExample('robotics/ComputeBSplineTrajectoryFor2DPlanarMotionExample')
% ----------------------------------------------------------------------
% Plot the results. Show the original waypoints, the computed polygon, and the interpolated B-spline.
figure;
hold all
plot(wptsx, wptsy, 'ob');
plot(xx1, yy1, '-k');
plot(cpts(1,:), cpts(2,:), '-sb');
plot(0.5*(cpts(1,1:end-1)+cpts(1,2:end)), 0.5*(cpts(2,1:end-1)+cpts(2,2:end)), '+b');
plot(q3(1,:), q3(2,:),'-*r');
legend('Original waypoints', 'Cubic spline', 'Computed control polygon', 'B-Spline');
% we get a quite similar curve for the spline, using pp instead of B-Spline
% for the visualization
% ----------------------------------------------------------------------
%% Interpretation of the LES to be solved to obtain the control polygon
% Visualization of a geometric method
% using lines and mid point construction
% solution is probably documentated somewhere, but not easy to be found
% ----------------------------------------------------------------------
% mid points of control polygon corners
midpts=0.5*(cpts(:,1:end-1)+cpts(:,2:end));
% mid points of mid points of control polygon corners 
midmidpts=0.5*(midpts(:,1:end-1)+midpts(:,2:end));

v1=midmidpts(:, 1:6)-cpts(:, 2:7);
v2=midmidpts-wpts;
d1=v1.^2;
d2=v2.^2;

l1=sqrt(d1(1,:)+d1(2,:));
l2=sqrt(d2(1,:)+d2(2,:));

% the ratio of lengthes is 1:3
% referencing to the midmid points.
l1./l2
% At least this is true for all interior points (unaffected by bc)

% the ratio of lengthes is 1:2
% referencing to the way points.
% This confirms the interpretation derived in the previous block
% At least this is true for all interior points (unaffected by bc)
figure; hold on;
    plot(wptsx, wptsy, 'ok');
    plot(cpts(1,:), cpts(2,:), '-sb');
    plot(q3(1,:), q3(2,:),'-*r');
    plot(midpts(1,:),midpts(2,:), '--+b');
 
    for k =1:6

        line([midmidpts(1,k),cpts(1,k+1)],[midmidpts(2,k),cpts(2,k+1)],'Color','m');
        line([midmidpts(1,k),wptsx(k)],[midmidpts(2,k),wptsy(k)],'Color','g');
    end
    legend('Original waypoints',  'Computed control polygon', 'B-Spline (low-res)', 'mid pts','2xVector1','Vector1');
% we get a quite similar curve for the spline, using pp instead of B-Spline

%% 
% Other Ex. 1: Draw spline within a given polygon
% Code reprint from The MathWorks Inc. help center
% ----------------------------------------------------------------------
points = [0 0 1 1 0 -1 -1  0   0 ;
          0 0 0 1 2 1   0 -2 -4]; 
 values = spcrv(points,3); 
 
figure;hold on;
    plot(points(1,:),points(2,:),'-o') ;
    plot(values(1,:),values(2,:),'-sk'), hold off;

%%
%  Other Ex. 2: Draw circular splines through points
% Code reprint from The MathWorks Inc. help center
% ----------------------------------------------------------------------
c1 = rscvn([-1 1 -1;0 0 0],[1 1;0 0]);
c = rscvn([-1 1;0 0],[1 1;0 0]);
% 180degree acrs 
% first from 0..pi
% second from pi to 2*pi=0
% slope at endpt is always x'=1, y'=0
figure;
    fnplt(c);
    axis([-1.05 1.05 -1.05 1.05]), axis equal;
    [form, order, breaks] = fnbrk(c,'f','o','b');

%%
% Other Ex. 3: Conversion to other types fails for full circle approx.
% Tested for 1 pt pp-spline and 2 pt B-Spline
% Seems doeable for half a circle, but correct bc's need to be found first
% Code reprint from The MathWorks Inc. help center
% ----------------------------------------------------------------------
pp2 = fn2fm(c,'B-',0.1); %two  pieces
pp1 = fn2fm(c,'pp',0.1); %just one piece, tribbled to 3 at bc
% interesting: we get a third dimension added
pp2.coefs;
pp1.coefs;

pp2.knots;
pp1.breaks;

% SURPRISE: something is added in the 3rd dimension !
% figure; %in 3D
% fnplt(pp1);
% %fnplt(pp2);
% axis([-1.05 1.05 -1.05 1.05])
 
yy = fnval(c, 0:0.1:2);
yy1 = fnval(pp1, 0:0.1:2);
yy2 = fnval(pp2, 0:0.1:2);

% The bc constraints are not converted
% We miss most of the circle
% plot the arc that is not a circle
figure; hold on;
    plot(yy(1,:),yy(2,:),'-k');
    plot(yy1(1,:),yy1(2,:),'-b');
    plot(yy2(1,:),yy2(2,:),':r');
    axis([-1.05 1.05 -1.05 1.05]), axis equal, grid on;

% the first dimension has been used as parameter i.e.
% yy(1,:) is linearly spaced
% yy(21,:) is quadratic, like y=1-0.5x^2 ~sqrt(1-x^2)
% what is the 2nd dimension? the y coords.

figure; hold on;
    plot(yy1(2,:),'-b');
    plot(yy2(2,:),':r');
    x=(0:0.05:1)*10+11;
    plot(x,0.5*(1-0.01*(x-11).^2),'-g');

%% Other Ex. 4: Some nice computer art
% Bronze Triskele Medallion in the Ulster Museum in Belfast
% Code reprint from The MathWorks Inc. help center 
% ----------------------------------------------------------------------
pp =[zeros(1,7); 5.4, 3, 6.9, 2.75, 2.5, .5, 5];
alpha = 2*pi/3; ca = cos(alpha); sa = sin(alpha); c = [ca sa;-sa ca];
d = [0 0 .05 -.05;1 -1 .98 .98]; d = [d c*d];
yin = rscvn([pp(:,[7,1:3]),c*pp(:,3:4),pp(:,3)], d(:,[1 2 1 4 7 5 1]));
figure;
    fnplt(yin), hold on, fnplt(fncmb(yin,c)), fnplt(fncmb(yin,c'))
    yang = rscvn([pp(:,6),-pp(:,6),pp(:,5),c*pp(:,4)],[d(:,[2 1 1]),c(:,2)]);
    fnplt(yang), fnplt(fncmb(yang,c)), fnplt(fncmb(yang,c'))
    axis([-7.2 7.2 -7.2 7.2]), axis equal, axis off, hold off