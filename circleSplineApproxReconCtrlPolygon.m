% Spline approximation of circles in GNU octave
% Reconstruction of the control polygon 
% ---------------------------------------------------------------------------
% Fit a piecewise cubic spline (of type pp) with breaks
% pp = splinefit (x, y, breaks)
% 
% From such an approximation of an circle, we want to reconstruct
% the outisde control polygon circumferencing the ellipse.
% 
% Using the derived control polygon, we can define a similar but
% B-Spline based approximation of the initial circle.
%
% The reconstruction of the outisde control polygon
% serves the purpose of converting a simple pp spline to a B-Spline

% Define an circle
% See splineExamples.m for details

% ---------------------------------------------------------------------------
% Definition of mathematicaly perfect ellipse:
% (x/a)^2+(y/b)^2=1
% ---------------------------------------------------------------------------
close all;
clear all;

% Abbreviations
% bc = boundary conditions

% Parameters
R=1;             % Squeeze a circle
skewness=1;    % into a ellipse of height=width*skewness
noPts = 5;      % number of sampling points (pick odd)

noPts=ceil(noPts/2)*2-1;

% Derived quantites
a=1*R;          % horizontal radius
b=skewness*R;   % vertical radius
R=[a,b];        % put everything into vectors
P0=[0,0];       % ellipse centre coordinates

% circumference of the ellipse
% for a mathematicaly perfect ellipse
circumference=2*sqrt(0.5*a^2+0.5*b^2)*pi;
% any straight line approximation will have a shorter circumference.


% ELLIPSE APPROXIMATION
% ---------------------------------------------------------------------------

% Griding to equidistant polar angles (low resolution sampling)
t1=linspace(-pi/2,3/2*pi,noPts);
V1=[cos(t1); sin(t1)]';
sampledEllipse=P0+R.*V1;
triangle1=diff(sampledEllipse,1);
d1 = sqrt(triangle1(:,1).^2+triangle1(:,2).^2);
p1=[0; cumsum(d1)];

% Derived quantites
curvature1 = (a*b)./((a*sin(t1)').^2+(b*cos(t1)').^2).^(3/2);
ROC=1./curvature1;% Radius of local curvature
maxL=5;           % max segment length for flat segments
maxAngle=15;      % increment in degree for circular segments
bowL=ROC*pi/180*maxAngle;
diffGrid=min(maxL, bowL); % avoid infinite ROC, cutoff at maxL

% Griding to equidistant polar angles (high resolution sampling)
t2=linspace(-pi/2,3/2*pi,101);
V2=[cos(t2); sin(t2)]';
highresEllipse=P0+R.*V2;
curvature2 = (a*b)./((a*sin(t2)').^2+(b*cos(t2)').^2).^(3/2);
triangle2=diff(highresEllipse,1);
d2 = sqrt(triangle2(:,1).^2+triangle2(:,2).^2);
p2=[0; cumsum(d2)];

% Griding to equidistant y-separations
% this will also produce more nodes in high curvture regions
t3=linspace(-1,1,(noPts-1)/2+1);
V3=[sqrt(1-t3(1:end-1).^2) -fliplr(sqrt(1-t3.^2)); t3(1:end-1) fliplr(t3)]';

equiYellipse=P0+R.*V3;
triangle3=diff(equiYellipse,1);
d3 = sqrt(triangle3(:,1).^2+triangle3(:,2).^2);
p3=[0; cumsum(d3)];

samplingMethod=1;
  switch (samplingMethod)
   case 1
    fromEllipse=sampledEllipse;
    fromGrid=p1;
   case 2
    fromEllipse=equiYellipse;
    fromGrid=p3;  
  end

  % May need to skip points around [-1,0] for error calculation
% Proper bc on 1st and 2nd derivatives would require 
% +2 data points at start and end of chain or use of periodic bc
% Simply skip first and last quarter of data to get rid of boundary impact
M=length(fromEllipse);
M4=(M-1)/4+1; idx=M4:M-M4-1;

% Local side length of triangle segments
triangle4=diff(fromEllipse,1);
% length of line segments - non-equal segment lengths
d4 = sqrt(triangle4(:,1).^2+triangle4(:,2).^2);
printf('Median segment length for %1.f pt approx. : %.3f\n', noPts, median(d4(idx)));
printf('Std dev of segment length for %1.f pt approx. : %.3f\n', noPts, std(d4(idx)));

% Assuming adequately sampled x,y data
% we can reparameterize into segments of approximately equal length
% using cubic spline interpolation
% parameter is the position along the circumference of the ellipse
% circumference=2*sqrt(0.5*a^2+0.5*b^2)*pi;
tm = cumsum([0; d4]);
approxCircumference=tm(end); % sum of straight line segments

% Define equidistant grid for interpolation 
% reconstruct high-resolution from downsampled optiEllipse
tt = linspace(0,approxCircumference,101);

% Interpolation of x and y values is done independently
% For height<width i.e. skewness<1 the flat segments at [0,+1] and [0,-1]
% will be nicely interpolated
% Large errors expected for the high curvature segments at [+1,0] and [-1,0]
% use periodic bc for full ellipse
px= splinefit(fromGrid,fromEllipse(:,1), fromGrid, "order", 3, "periodic", true); 
py= splinefit(fromGrid,fromEllipse(:,2), fromGrid, "order", 3, "periodic", true); 
splineEllipse(:,1) =ppval(px,tt);
splineEllipse(:,2) =ppval(py,tt);
% for linear use
%splineEllipse=interp1(fromGrid,fromEllipse,tt,'linear');
%test and compared to interp1 using "linear", "spline", "pchip"

triangleNormalized=(splineEllipse./R).^2;
lengthNormalized=sqrt(triangleNormalized(:,1).^2+triangleNormalized(:,2).^2);
% this is the error distance to the perfect ellipse
%err=sqrt((1-lengthNormalized.^2));
err=1-lengthNormalized;
% Need to skip points around [-1,0] for error calculation
% Proper bc on 1st and 2nd derivatives would require 
% +2 data points at start and end of chain
MM=length(splineEllipse);
MM4=(MM-1)/4+1; jdx=MM4:MM-MM4-1;
rms=sqrt(sum(err(jdx).^2))/length(tt);
rms=(sum(err(jdx)))/length(tt);
% this is not the real rms error as usual
% it is just the missmatch with the implicit ellipse definition
% 1/N*sum(xi/a)^2+(yi/b)^2=1-rms(of average single point)


% Report analysis results
triangle=diff(splineEllipse,1);
d = sqrt(triangle(:,1).^2+triangle(:,2).^2);

printf('Median segment length for reconstructed 101 pt interpolation: %.3f\n', median(d(jdx)));
printf('Std dev of segment length reconstructed 101 pt interpolation:: %.3f\n', std(d(jdx)));
printf('RMS error for reconstructed interplation: %.3f\n', rms);
printf('Max. deviation from perfect ellipse is: %.3f\n', max(err(jdx)));

% ---------------------------------------------------------------------------
% Report figures for illustration
plotit=false;
if (plotit==true)
figsize=1600;
fig1=figure(1, 'position',figsize*[0,0,1,1]+[500,200,50,0]);

axis([-1 1 -1 1]);
axis equal
grid on
hold on
plot(highresEllipse(:,1),highresEllipse(:,2) ,'k');
%plot(optiEllipse(:,1),optiEllipse(:,2) ,'-ob');
plot(sampledEllipse(:,1),sampledEllipse(:,2) ,'-sg');
plot(splineEllipse(:,1),splineEllipse(:,2) ,'-m');
plot(equiYellipse(:,1),equiYellipse(:,2) ,'-c');
axis equal


% Attention: x and y values are interpolated independently from each other
% This is obviously not the best interpolation methode
% plot interpolation of x and y-values
figsize=1600;
fig2=figure(2, 'position',figsize*[0,0,1,1]+[500,200,50,0]+figsize*[1,0,0,0]);
hold on;
plot(p1/p1(end),sampledEllipse(:,1),'-sg');
plot(p2/p2(end),highresEllipse(:,1),'-k');
plot(tt/tt(end),splineEllipse(:,1),'-om');
%plot(optiGrid/optiGrid(end),optiEllipse(:,1),'-ob');
plot(p3/p3(end),equiYellipse(:,1) ,'-c');
plot(p1/p1(end),sampledEllipse(:,2),'-sg');
plot(p2/p2(end),highresEllipse(:,2),'-k');
plot(tt/tt(end),splineEllipse(:,2),'-om');
%plot(optiGrid/optiGrid(end),optiEllipse(:,2),'-ob');
plot(p3/p3(end),equiYellipse(:,2) ,'-c');
end
% ---------------------------------------------------------------------------
% END ELLIPSE APPROXIMATION

% ---------------------------------------------------------------------------
% RECONSTRUCT THE CTRL POLYGON
% See splineExamples.m for more details.

% ----------------------------------------------------------------------
% Create waypoints to interpolate with a B-Spline.

%wptsx = [fromEllipse(end-1,1) fromEllipse(:,1)' fromEllipse(2,1)];
%wptsy = [fromEllipse(end-1,2) fromEllipse(:,2)' fromEllipse(2,2)];

%wptsx = [fromEllipse(1,1) fromEllipse(:,1)' fromEllipse(end,1)];
%wptsy = [fromEllipse(1,2) fromEllipse(:,2)' fromEllipse(end,2)];


%t0=acos(0.7);
%V0=[cos(t1); sin(t1)]'
%P0=P0+R.*V1


%wpts=1.0*[0.6667    0.6667    0.0000   -1.0000   -0.0000    0.6667    0.6667;
%          0.0000    0.3491    0.9933    0.0000   -0.9933   -0.3491   -0.0000];

% Simple but not best approach to keep number of knots low
% We do pad the knot vector with a second 360 degree period
wpts =[    0  1  0  -1  0  1 0  ;   %x
          -1  0  1  0  -1  0 1   ]; %y
       
% Put in vector to solve later a linear equation system LES.
wptsx=wpts(1,:);
wptsx=wpts(2,:);
L = length(wpts) - 1;

% Form matrices used to compute INTERIOR points of control polygon
% Interior points are just the normal pts excluding the start and end pt
r = zeros(L+1, size(wpts,1));
A = eye(L+1);
for i= 2:L
    A(i,(i-1):(i+1)) = [1 4 1];  %1+4+1 = 6
    r(i,:) = 6*wpts(:,i)';       %tranposed vector wpts
end

%%
% If you know how to specifiy boundary conditions (bc's)
%   this is the section to show your brilliance.
%   Contribute with review comments - post your own methods.
%   Uploads and share your knowledge.
%
% We do not want to add periodic bc's
%A(1,1:2) = [4 1];
%A(1, L+1)  = 1; 
%A(L+1,L:L+1) = [1 4];
%A(L+1,1)  = 1; 

% We did pad the knot vector with a second 360 degree period
% What boundary values to patch?
% Figure it out if you can, post and share it.
A(2,1:3) = [3/2 7/2 1]; 
A(L,(L-1):(L+1)) = [1 7/2 3/2]; 
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

v1=midmidpts(:, 1:end)-cpts(:, 2:end-1);
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

% Plot the results. Show the original waypoints, the computed polygon, and the interpolated B-spline.
% Report figures for illustration
figsize=1600;
fig1=figure(1, 'position',figsize*[0,0,1,1]+[500,200,50,0]);
hold on;
axis([-1 1 -1 1]*2);
axis equal
grid on
hold on


plot(fromEllipse(:,1), fromEllipse(:,2), '-b');
plot(highresEllipse(:,1), highresEllipse(:,2), '-k');
skip=0;
plot(cpts(1,skip+1:length(midmidpts)-skip), cpts(2,skip+1:length(midmidpts)-skip), '-sr');

plot(midpts(1,skip+1:length(midmidpts)-skip),midpts(2,skip+1:length(midmidpts)-skip), '--+r');

    for k =skip+1:length(midmidpts)-skip

        line([midmidpts(1,k),cpts(1,k+1)],[midmidpts(2,k),cpts(2,k+1)],'Color','m');
        line([midmidpts(1,k),wptsx(k+1)],[midmidpts(2,k),wptsy(k+1)],'Color','g');
    end
%plot(q1(1,:), q1(2,:),'-r');

 legend('Original waypoints', 'pp Spline (high-res)','Computed control polygon','mid pts','2xVector1','Vector1');

% ---------------------------------------------------------------------------
% END RECONSTRUCT THE CTRL POLYGON