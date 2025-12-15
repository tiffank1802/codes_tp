function [u,v,a] = Newmark (M, C, K, f, u0, v0, a0, dt, beta, gamma, varargin)
% This function solves the dynamic system with a Newmark approach
% Ma + Cv + Ku = f

% input : f        : mulit-row matrix of the loading
%         M        : mass matrix
%         C        : Damping matrix
%         K        : stiffness matrix (with its Lagrange multipliers)
%         u0       : displacement at t = 0
%         v0       : velocity at t = 0
%         a0       : accelerations at t = 0 /!\ should be determined via Ma = f-Cv-Ku /!\
%         dt       : time discretization parameter (constant)
%         beta     : parameter of the method (centered accel: 1/4)
%         gamma    : parameter of the method (centered accel: 1/2)
 
% output : u    : multi-row matrix of the solution
%          v    : multi-row matrix of the velocity
%          a    : multi-row matrix of the acceletation

 ndof  = size(f,1);
 ntime = size(f,2);

 u = zeros( ndof, ntime );
 v = zeros( ndof, ntime );
 a = zeros( ndof, ntime );
 
 up = u0; vp = v0; ap = a0;
 %ap = M\(f(:,1)-C*vp-K*up); 

 fatK = K + 1/(beta*dt^2)*M + gamma/(beta*dt)*C;
 
% [Ll, Uu, Pp] = lu (fatK); % Speed things up
% [Lb,Ub,Pb,Qb] = lu (fatK,"vector");
 
 invert = 0;
 if numel(varargin) > 0
    fatK1 = varargin{1};
    invert = 1;
 elseif 2*ntime > ndof % It's preferable to invert fatK
    fatK1 = inv(fatK);
    invert = 1;
 else
    [Ll,Uu,Pp,Qq] = lu (fatK,"vector"); % TODO: see for a choelesky instead
 end
 
 u(:,1) = u0; v(:,1) = v0; a(:,1) = a0;
 
 for i=2:ntime
    b = ( f(:,i) + ...
          C*( gamma/(beta*dt)*up + (gamma/beta-1)*vp + dt/2*(gamma/beta-1)*ap ) + ...
          M*( 1/(beta*dt^2)*up + 1/(beta*dt)*vp + (1/(2*beta)-1)*ap ) );
          
    if invert
       u(:,i) = fatK1*b;
    else
       u(Qq,i) = Uu \ ( Ll \ (b(Pp)) );
    end
               
    a(:,i) = 1/(2*beta*dt^2/2) * ( u(:,i) - up - dt*vp - dt^2/2*(1-2*beta)*ap );
    v(:,i) = vp + dt*( (1-gamma)*ap + gamma*a(:,i) );
    
    up = u(:,i); vp = v(:,i); ap = a(:,i);
 end

end
