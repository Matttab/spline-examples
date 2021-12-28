% Spline approximation of ellipses in GNU octave
% How to choose appropriate nodes on the circumference
% ---------------------------------------------------------------------------
% Fit a piecewise cubic spline (of type pp) with breaks
% pp = splinefit (x, y, breaks)
% ---------------------------------------------------------------------------
% Spline approximation of a closed ellipse with a total of 9 nodes
% is already usable on common CNC xy milling machines. The approximation
% is better than for stepwise linear interpolation, but requires support
% of spline interpolation (G5, G5.1 or G5.2 commands).
%
%   G5 creates a cubic B-spline
%   G5.1 creates a quadratic B-Spline
%   G5.2 G5.3 NURBS Block rational B-Spline
%
% Spline interpolation yields good results for 9 or more nodes. 
% Definitions: Nodes = knots = waypoints of tool 
% This is true for various potential positions of the nodes on the
% circumference of a mathematicaly perfect ellipse.
% ---------------------------------------------------------------------------
% Tested are an equidistant distribution in polar angle space and a
% Method** novel and uncommon, with purpose to get a higher 
%          node density around regions of high curvature. **
% Also equidistant distribution along the short ellipse major axis is tested.
% This will also increase the density of nodes around high curvature parts.
% ---------------------------------------------------------------------------
% 
% For the minimum of 9 nodes, the newly proposed method**
% is slightly superior (on cost of some computation efforts).
% For 11 nodes and more, equidistant distribution in polar angle space is
% just equally fine. Most CAD, CAM or post-processers for CNC machining will
% simply transform any type of spline to short straight segments prior
% to sending it to the G-Code interpreter.
%
% A more simple approach than ** would slight shift the nodes
% towards the high curvature regions by applying a distortion to
% the equidistant polar grid. The optimum position of the 2 intermediates nodes
% per quadrant can be found numerically. A list of these positions 
% for ellipses with varying skewness would suit for most practial purposes.
%
% G-Code G1 straight line interpolation needs more nodes to get similarly
% low deviations from an perfect ellipse. A total of 21 nodes is likely
% sufficient for the most practically relevant CNC tasks. About 71 nodes would
% yield similarly good approximations as a spline approach with 13 nodes.

% The code of the published source
% https://www.mathworks.com/matlabcentral/answers/86615-how-to-plot-an-ellipse#answer_96122
% has been used as starting base.
%
% Newly introduced was a vector format, error calculation, and most interestingly,
%    a novel griding methode**
%    respecting he varying curvature along the circumference of the ellipse.
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
skewness=0.3;    % into a ellipse of height=width*skewness
noPts = 13;      % number of sampling points (pick odd)

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

% ---------------------------------------------------------------------------
% Method** content
% generate a grid with a local spacing of diffGrid around nodes of grid fromGrid
%      use some approximation to achieve that
%      here we tak a -cosine function for the spacings
%      the integral optiGrid4 is a -sinus function that can be directly generated
%      diff(optiGrid4) is somewhat similar to the required variable spacing
%      having same min and max spacing, a smooth crossover, and equidistant
%      segment length around the x and y axis interceptions (angle = integer * pi/2)
A=max(diffGrid)-min(diffGrid);
const=(max(diffGrid)+min(diffGrid))/2;

% interpolmate half of the ellipse
N4 = (noPts-1)/4;   % no of segments in one quadrant
N=2*(2*N4+1)-1;     % no of segments for full ellipse
ttt4=[0:1/N4:2]*pi;

% define the integral directly as sinus function
% the derivates will be a cosine as targeted
% it may need to fix the phase and sign to get the targeted spacing
optiGrid4=0.5*A*N4/pi*sin(ttt4)+N4/pi*const*ttt4;

% interpolmate full of the ellipse
ttt=[-1/2:1/(N-1)*2:3/2]'*pi;
optiGrid1=[optiGrid4(1:end-1) optiGrid4(1:end)+optiGrid4(end)-optiGrid4(1)];

% scale to correct circumference
% optiGrid=optiGrid1/optiGrid1(end)*circumference;
% use optiGrid to grid the circumference instead of equidistant gridding
% the optiGrid has more dense sampling at high curvature segments
% optiGrid1 is still in mm units along the circumference
% figure; hold on;
% plot(t1(1:end-1), diff(p3),'b')
% plot(ttt(1:end-1), diff(optiGrid1),'r')

% Improved griding respecting the variable curvature
optiEllipse = interp1(p2,highresEllipse,optiGrid1/optiGrid1(end)*p2(end),"linear");
% you can test and compare 'linear', 'spline', 'pchip'
optiTriangle=diff(optiEllipse,1);
optiD = sqrt(optiTriangle(:,1).^2+optiTriangle(:,2).^2);
optiP=[0; cumsum(optiD)];
optiGrid=optiGrid1/optiGrid1(end)*optiP(end);

% End Method** content
% ---------------------------------------------------------------------------

% compare reconstruction from
%   p1, sampledEllipse vs
%   optiGrid, optiEllipse
% select wich method to compare, report rms error, and plot
% Reported is not the real rms error as usual
% it is just the missmatch with the implicit ellipse definition
%   1/N*sum(xi/a)^2+(yi/b)^2=1-rms(of average single point)
% serving fine as figure of merit.

samplingMethod=2;
  switch (samplingMethod)
   case 1
    fromEllipse=sampledEllipse;
    fromGrid=p1;
   case 2
    fromEllipse=optiEllipse;
    fromGrid=optiGrid;
   case 3
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
printf('Median segment length for low-res sampling: %.3f\n', median(d4(idx)));
printf('Std dev of segment length low-res sampling: %.3f\n', std(d4(idx)));

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
err=sqrt((lengthNormalized-1).^2);
% Substitute with err=1-lengthNormalized and recalc result.

% Need to skip points around [-1,0] for error calculation
% Proper bc on 1st and 2nd derivatives would require 
% +2 data points at start and end of chain
MM=length(splineEllipse);
MM4=(MM-1)/4+1; jdx=MM4:MM-MM4-1;
rms=sqrt(sum(err(jdx).^2))/length(tt);
% Substitute with rms=(sum(err(jdx)))/length(tt); and recalc result.

% this is not the real rms error as usual
% it is just the missmatch with the implicit ellipse definition
% 1/N*sum(xi/a)^2+(yi/b)^2=1rms(of average single point)


% APPROXIMATION RESULT SUMMARY
% ---------------------------------------------------------------------------
% RMS error example for skewness=0.3, 101 pts generated from 9 pts straight-seg trace

 % Result for sampling of equidistant polar angles
  %   Median segment length: 0.382+/-0.130 (by definition equidistant in polar angle)
  %   Spline interpolat166ion: 0.011, Max. err: 0.296
  % Result for sampling of optimized sampling according to inverse curvature
  %   Median segment length: 0.382+/-0.244 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Spline interpolation: 0.009, Max. err: 0.254

% RMS error example for skewness=0.3, 101 pts generated from 13 pts straight-seg trace
  
  % Result for sampling of equidistant polar angles
  %   Median segment length: 0.382+/-0.130 (by definition equidistant in polar angle)
  %   Spline interpolat166ion: 0.012, Max. err: 0.267
  % Result for sampling of optimized sampling according to inverse curvature
  %   Median segment length: 0.382+/-0.244 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Spline interpolation: 0.0
  % Result for sampling of equidistant y coordinates
  %   Median segment length: 0.221+/-0.313 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Spline interpolation: 0.013, Max. err: 0.297
  
  % RMS error example for skewness=0.3, 101 pts generated from 21 pts straight-seg trace
  % Result for sampling of equidistant polar angles
  %   Median segment length: 0.231+/-0.077 (by definition equidistant in polar angle)
  %   Spline interpolation: 0.012, Max. err: 0.295
  % Result for sampling of optimized sampling according to inverse curvature
  %   Median segment length: 0.218+/-0.148 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Spline interpolation: 0.012, Max. err: 0.295
  % Result for sampling of equidistant y coordinates
  %   Median segment length: 0.221+/-0.313 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Spline interpolation: 0.012, Max. err: 0.292
  
% Almost equivalent result!
% The ellipse can be reconstructed with splines equivally well from
% a low-reolution sample data
%    - based on equidistant parameterization of polar angle
%    - based on non-equidistant parameterization according to inverse of local curvature

% Let's compare same grids for straight linear interpolation (e.g. if no spline interpolation is available) 

  % RMS error example for skewness=0.3, 101 pts generated from 13 pts straight-seg trace
  % Result for sampling of equidistant polar angles
  %   Median segment length: 0.382+/-0.129 (by definition equidistant in polar angle)
  %   Spline interpolation: 0.015, Max. err: 0.340
  % Result for sampling of optimized sampling according to inverse curvature
  %   Median segment length: 0.382+/-0.244 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Spline interpolation: 0.016, Max. err: 0.337
  % Result for sampling of equidistant y coordinates
  %   Median segment length: 0.221+/-0.313 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Spline interpolation: 0.017, Max. err: 0.357
  
  % RMS error example for skewness=0.3, 101 pts generated from 21 pts straight-seg trace
  % Result for sampling of equidistant polar angles
  %   Median segment length: 0.231+/-0.077 (by definition equidistant in polar angle)
  %   Linear interpolation: 0.013, Max. err: 0.310
  % Result for sampling of optimized sampling according to inverse curvature
  %   Median segment length: 0.218+/-0.148 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Linear interpolation: 0.014 Max. err: 0.313
  % Result for sampling of equidistant y coordinates
  %   Median segment length: 0.221+/-0.313 (by purpose non-equidistant, distance scaling with inverse curvature)
  %   Spline interpolation: 0.014, Max. err: 0.307

% Almost equivalent result!
% The rms error is similar for both griding methods.

% END APPROXIMATION RESULT SUMMARY
% ---------------------------------------------------------------------------


% Report analysis results
triangle=diff(splineEllipse,1);
d = sqrt(triangle(:,1).^2+triangle(:,2).^2);

printf('Median segment length for reconstructed 101 pt interpolation: %.3f\n', median(d(jdx)));
printf('Std dev of segment length reconstructed 101 pt interpolation:: %.3f\n', std(d(jdx)));
printf('RMS error for reconstructed interplation: %.3f\n', rms);
printf('Max. deviation from perfect ellipse is: %.3f\n', max(err(jdx)));

% ---------------------------------------------------------------------------
% Report figures for illustration
figsize=1600;
fig1=figure(1, 'position',figsize*[0,0,1,1]+[500,200,50,0]);

axis([-1 1 -0.5 1.5]);
axis equal
grid on
hold on
plot(highresEllipse(:,1),highresEllipse(:,2) ,'k');
plot(optiEllipse(:,1),optiEllipse(:,2) ,'-ob');

plot(splineEllipse(:,1),splineEllipse(:,2) ,'-m');
plot(sampledEllipse(:,1),sampledEllipse(:,2) ,':sg');

plot(equiYellipse(:,1),equiYellipse(:,2) ,'--+c');
leg=legend({'high-res polar angle interpolation','nodes shifted to high curvatures', ...
'splinefit to low-res equidistant polar angle','low-res equidistant polar angle interpolation', ...
'low-res equidistant y separation'},'location','northwest')
   set (leg, "fontsize", 14);
   set(leg,'position',[0.25 0.5 0.55 0.3])
%axis equal


% Attention: x and y values are interpolated independently from each other
% This is obviously not the best interpolation methode
% plot interpolation of x and y-values
figsize=1600;
fig2=figure(2, 'position',figsize*[0,0,1,1]+[500,200,50,0]+figsize*[1,0,0,0]);
hold on;
plot(p2/p2(end),highresEllipse(:,1),'-k');
plot(optiGrid/optiGrid(end),optiEllipse(:,1),'-ob')
plot(tt/tt(end),splineEllipse(:,1),'-m');
plot(p1/p1(end),sampledEllipse(:,1),':sg');
plot(p3/p3(end),equiYellipse(:,1) ,'--+c');


plot(p2/p2(end),highresEllipse(:,2),'-k');
plot(optiGrid/optiGrid(end),optiEllipse(:,2),'-ob');
plot(tt/tt(end),splineEllipse(:,2),'-m');
plot(p1/p1(end),sampledEllipse(:,2),':sg');
plot(p3/p3(end),equiYellipse(:,2) ,'--+c');
grid on
leg=legend({'high-res polar angle interpolation','nodes shifted to high curvatures', ...
'splinefit to low-res equidistant polar angle','low-res equidistant polar angle interpolation', ...
'low-res equidistant y separation'},'location','northwest')
   set (leg, "fontsize", 12);
   set(leg,'position',[0.5 0.65 0.45 0.3])
xlabel('Normalized distance along circumference. Starting at (x,y)=[0,-0.3]');
ylabel('x (sin-like) and y (cos-like) coordinates');
% ---------------------------------------------------------------------------
% END FILE