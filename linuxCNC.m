classdef linuxCNC < handle
% Some routines taken from linuxCNC source code to interpolate splines
% Code converted from C++ to ocatve (may code translation errors now)
% 
% Math reference:
% https:%en.wikipedia.org/wiki/B%C3%A9zier_curve


properties
  metric=true;
   _pos_x;
   _pos_y;
   _pos_z;
   _pos_a;
   _pos_b;
   _pos_c;
   _pos_u;
   _pos_v;
   _pos_w;
end %properties

methods (Static=true)
  function bersteinCoeff=kanonical2bernstein(kanonicalCoeff)
  M=[1 -3 3 -1; %1-3t+3t^2- t^3
     0 3 -6 3;  %  3t-6t^2+3t^3
     0 0 3 -3;  %     3t^2-3t^3
     0 0 0  1]; %          3t^3
  bersteinCoeff=((M*kanonicalCoeff')');
  end
  function kanonicalCoeff=bernstein2kanonical(bersteinCoeff)
  Mi=[3 3 3 3;
      0 1 2 3;
      0 0 1 3;
      0 0 0 3]/3;
  kanonicalCoeff=((Mi*bersteinCoeff')');
  end
end %end static methods

methods 
  
% [x,y]=SPLINE_FEED2( obj, x1,  y1,  x2,  y2) 
% End point = [x2,y2]
% Start point = [obj._pos_x,obj._pos_y]
% Tanget slope is (y1-y0)/(x1-x0) at start
% Tanget slope is (y2-y1)/(x2-x1) at end
function [x,y]=SPLINE_FEED2( obj, x1,  y1,  x2,  y2) 
     x0 = obj._pos_x; %start point
     y0 = obj._pos_y; %start point
    % Tanget slope is (y1-y0)/(x1-x0) at start
    % Tanget slope is (y2-y1)/(x2-x1) at end
    fprintf(stderr, "SPLINE_FEED(quadratic): %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
          x0,y0,x1,y1,x2,y2);

    for i=0:100 % do 100 straight line per segment  
       t = i / 100.;
       % quadratic Bezier base
       % Bk2 = (2 k)t^k*(1-t)^(2-k)
       B22 = t*t;        % dB22/dt=2*t=         [0   2]     for t=[0 1]
       B12 = 2*t*(1-t);  % dB12/dt=2(1-2*t)=    [2  -2]  for t=[0 1]
       B02 = (1-t)*(1-t);% dB02/dt=-2*(1-t)=    [-2 -0]  for t=[0 1]
       x(i+1,1) = x0*B02 + x1*B12 + x2*B22;%= 2*[(x1-x0) (x2-x1)]
       y(i+1,1) = y0*B02 + y1*B12 + y2*B22;%= 2*[(y1-y0) (y2-y1)]
       %dx/dy=(dx/dt)/(dy/dy)=(y1-y0)/(x1-x0) at t=0
       %dx/dy=(dx/dt)/(dy/dy)=(y2-y1)/(x2-x1) at t=1
      % EBo -- replace 12345 with *whatever* gives us the line_number
      % STRAIGHT_FEED(12345, x,y, _pos_z, _pos_a, _pos_b, _pos_c, _pos_u, _pos_v, _pos_w);
    end
end

% [x,y]=SPLINE_FEED2( obj, x1,  y1,  x2,  y2, x3, y3) 
% End point = [x3,y3]
% Start point = [obj._pos_x,obj._pos_y]
% Tanget slope is (y1-y0)/(x1-x0) at start
% Tanget slope is (y3-y2)/(x3-x2) at end
function [x,y]=SPLINE_FEED3( obj, x1,  y1,  x2,  y2,  x3,  y3) 
     x0 = obj._pos_x;
     y0 = obj._pos_y;
    % Tanget slope is (y1-y0)/(x1-x0) at start
    % Tanget slope is (y3-y2)/(x3-x2) at end
    fprintf(stderr, "SPLINE_FEED(cubic): %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
          x0,y0,x1,y1,x2,y2,x3,y3);
          
    for i=0:100 % do 100 straight line per segment  
      % cubic Bezier base for order n=3
      % bi,n = (n k)*t^i*(1-t)^(n-i)
      % start point interpolation property: start at t1(0)=1
      % end   point interpolation property: end   at t3(0)=1
      % Bk3 = (3 k)t^k*(1-t)^(3-k)
       t = i / 100;
       B33 = t*t*t;            % t3(0)=0[3x] and                 t3(1)=1[3] start point
       B23 = 3*t*t*(1-t);      % t3(0)=0[2x] and t3(1)=0[1x]                             t2(1/2)=3/8
       B13 = 3*t*(1-t)*(1-t);  % t3(0)=0[1x] and t3(1)=0[2x]                             t2(1/2)=3/8
       B03 = (1-t)*(1-t)*(1-t);%             and t3(1)=0[3x] and t3(0)=1[3] endpoint     t2(1/2)=3/8
       x(i+1,1) = x0*B03 + x1*B13 + x2*B23 + x3*B33;
       y(i+1,1) = y0*B03 + y1*B13 + y2*B23 + y3*B33;
       
       %B33/dt=3*t^2 =        [0    3]     for t=[0 1] B33'= 3*[B22]=3*t^2 
       %B13/dt=3t*(2-3*t)=    [0   -3]     for t=[0 1] B13'= 3*[B02-B12]=3t*(2-3*t)
       %B23/dt=3*(1-3t)(1-t)= [3    0]     for t=[0 1] B23'= 3*[B12-B22]=3*(1-t)(1-3*t)
       %B03/dt=-3(1-t)^2=     [-3   0]     for t=[0 1] B03'= 3*[-B02]=-3(1-t)^2
       %dx/dy=(dx/dt)/(dy/dy)=(y1-y0)/(x1-x0) at t=0
       %dx/dy=(dx/dt)/(dy/dy)=(y3-y2)/(x3-x2) at t=1
       
       
      % EBo -- replace 12345 with *whatever* gives us the line_number
      % STRAIGHT_FEED(12345, x,y, _pos_z, _pos_a, _pos_b, _pos_c, _pos_u, _pos_v, _pos_w);
    end
end
% 
% C=[ 1 2 3 4; 1 1.4 1.5 1.55]'
% pts= construct_spline(C)
function [x, y]=construct_spline(obj,input)
   C=input';
    %4 coeffs. for cubic spline a=t0,b=t1,c=t2,d=t3
    n = size(C,1)-1; % number of splines
    S=zeros(n,2,4);   % 2dim = x,y, 
    for dim=1:2  %x,y loop
        % define d:
          d(1)=3.0*(C(2,dim)-C(1,dim));
        for i=2:n
          d(i)=3.0*(C(i+1,dim)-C(i-1,dim));
        end
          d(n+1)=3.0*(C(n+1,dim)-C(n,dim));
          %  forward sweep: simplified tridiagonal solution:
          % ball a(i)=1 and all c(i)=1
        
            b(1)=2.0;
        for i=2:n+1
            w = 1.0/b(i-1);
            b(i) = 4.0 - w;
            d(i) -= (w*d(i-1));
        end
        % calculate solution vector x(i) = D(i):
        % (Wikipedia x() = Wolfram D())
        x(n+1)=d(n+1)/b(n+1);
        for i=n:-1:1
            x(i)=(d(i)-x(i+1))/b(i);
         end
        % calculate spline S(i)(dim) a, b, c and d:
        for i=1:n
            S(i,dim,1)=C(i,dim);
            S(i,dim,2)=x(i);
            S(i,dim,3)=3.0*(C(i+1,dim) -C(i,dim)) - 2.0*x(i) - x(i+1);
            S(i,dim,4)=2.0*(C(i,dim)   -C(i+1,dim)) + x(i) + x(i+1);
        end
    end
        for p=1:n  %spline points
    
          for k=0:11  %time 0-1 for each spline
            t=k*0.1; 
             
             % x,y = d*t^3 + c*t^2 + b*t + a;
            x(p,k+1) = ((S(p,1,4)*t +S(p,1,3))*t +S(p,1,2))*t +S(p,1,1);
            y(p,k+1) = ((S(p,2,4)*t +S(p,2,3))*t +S(p,2,2))*t +S(p,2,1);
            %std::cout<< "x:" << px << " y:" << py << std::endl;          
          end
        end

end
    

% template for NURBS_FEED
% Code not functional
function NURBS_FEED(obj, nurbs_control_points, k) 
     u = 0.0;
     n = length(nurbs_control_points) - 1;
     umax = n - k + 2;
     div = (n+1)*15;
     knot_vector = zeros(n,k);	
     P1 =[];
    while (u+umax/div < umax) 
        P1 = nurbs_point(u+umax/div, k, nurbs_control_points, knot_vector);
        % EBo -- replace 12345 with *whatever* gives us the line_number
        %STRAIGHT_FEED(12345, P1.X,P1.Y, _pos_z, _pos_a, _pos_b, _pos_c, _pos_u, _pos_v, _pos_w);
        u = u + umax/div;
    end 
    P1(1) = nurbs_control_points(1,n);
    P1(2) = nurbs_control_points(2,n);
    % EBo -- replace 12345 with *whatever* gives us the line_number
    % STRAIGHT_FEED(12345, P1.X,P1.Y, _pos_z, _pos_a, _pos_b, _pos_c, _pos_u, _pos_v, _pos_w);
    knot_vector=[];
end

end %end methods



end %class