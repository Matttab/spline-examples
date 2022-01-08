clear all; close all;
mySpline = linuxCNC;

mySpline.metric=true;

% the control points of the B-spline 
% are asymmetrically placed - SURPRISE !
% control points are asymmetric to get a symmetric shape

unitsphere = [1  0 -1  0  1;  %x
              0  1  0 -1  0]; %y
       
  
 unitcoefs = ...
 [    1.000    0.520  -0.520  -1.000   -1.000  -0.520  0.520  1.000;
      0.520    1.000   1.000   0.520   -0.520  -1.000 -1.000 -0.5200];




% Parameters
R=10;             % Squeeze a circle
skewness=0.5;    % into a ellipse of height=width*skewness
noPts = 101;      % number of sampling points (pick odd)

noPts=ceil(noPts/2)*2-1;

% Derived quantites
a=1*R;          % horizontal radius
b=skewness*R;   % vertical radius
R=[a,b];        % put everything into vectors
P0=[0,0];       % ellipse centre coordinates

% squeeze circle to ellipse
 pts=(R.*unitsphere')' ;
 coefs=(R.*unitcoefs')' ;


filename='splineEllipseTest.ngc';
fileID = fopen(filename,'w+'); % open a file and overwrite content

 
splinesegments=[];
for k=1:4% 1:length(ctrlPt)

 % copy start and end pts to be used as ctr pts
ctrl3Pt(:,1)=pts(:,k);    
ctrl3Pt(:,4)=pts(:,k+1);
 x0=pts(1,k);
 y0=pts(2,k);
 mySpline._pos_x=x0;
 mySpline._pos_y=y0;


  % spline feed
  % optimize location of ctrl pts 
  %ctrlPt=[0.551:0.001:0.553];   %use 0.552
  ctrlPt=0.552;
  % 2 ctr points per segment
  idx=2*(k-1)+1;
  ctrl3Pt(:,2)=coefs(:,idx); 
  ctrl3Pt(:,3)=coefs(:,idx+1);   
  ctrl3Pt
  [xc,yc,splinecode]=mySpline.SPLINE_FEED3(ctrl3Pt(1,2), ctrl3Pt(2,2),   ctrl3Pt(1,3),  ctrl3Pt(2,3),  ctrl3Pt(1,4), ctrl3Pt(2,4) ); 
  %d=diff(yc)./diff(xc);
  %di=diff(xc)./diff(yc); %inverse slope
  fdisp(fileID, splinecode);
  splinecode=[];
  splinesegments{k}=[xc,yc];

  end

fclose(fileID); 

  % model ideal curve
  xc=splinesegments{1}(:,1);
  yc=splinesegments{1}(:,2);
  no=length(xc);
  phi1=atan(yc./xc/skewness)';
  rho1=sqrt(yc.^2+xc.^2);
  V1=[cos(phi1); sin(phi1)]';
  idealCurve=P0+R.*V1;
  err=(yc-idealCurve(:,2)).^2+(xc-idealCurve(:,1)).^2;
  rms(k)=sqrt(sum(err)/no)

  
%[minval,minidx]=min(rms);
%ctrlPt(minidx)
%figure
%plot(ctrlPt, rms)

% circumference of the ellipse
% for a mathematicaly perfect ellipse
circumference=2*sqrt(0.5*a^2+0.5*b^2)*pi;
% any straight line approximation will have a shorter circumference.

% Griding to equidistant polar angles (low resolution sampling)
t2=linspace(0,2*pi,noPts);
V2=[cos(t2); sin(t2)]';
sampledEllipse=P0+R.*V2;
triangle1=diff(sampledEllipse,1);
d2 = sqrt(triangle1(:,1).^2+triangle1(:,2).^2);
p2=[0; cumsum(d2)];


figure; hold on;
plot(sampledEllipse(:,1), sampledEllipse(:,2),'--c');
plot(pts(1,:), pts(2,:),':+');
plot(coefs(1,:), coefs(2,:),'--o');
axis([-1 1 -1 1]*max(R)*1.1);

plot(idealCurve(:,1),idealCurve(:,2),'-b');
for k=1:4
plot(splinesegments{k}(:,1),splinesegments{k}(:,2),'-r');
end
axis square;
grid on;

