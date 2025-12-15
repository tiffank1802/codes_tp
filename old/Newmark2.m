function [u,v,a] = Newmark2 ( M, C, K, f, u0, v0, dt, beta, gamma )
% This function solves the dynamic system with a Newmark approach (incremental improvement on more basic Newmark.m)
% Ma + Cv + Ku = f

% input : f        : multi-row matrix of the loading
%         M        : mass matrix
%         C        : Damping matrix
%         K        : stiffness matrix (with its Lagrange multipliers)
%         u0       : displacement at t = 0
%         v0       : velocity at t = 0
%         dt       : time discretization parameter (constant)
%         beta     : parameter of the method (centered accel: 1/4)
%         gamma    : parameter of the method (centered accel: 1/2)
%         fatK1    : invert of fatK (optionnal)
%         ineq     : A structure for imposing an inequality Cu<d via Uzawa
 
% output : u    : multi-row matrix of the solution
%          v    : multi-row matrix of the velocity
%          a    : multi-row matrix of the acceletation

 fatK1 = false; ineq = false; % De-activate funky functions. TODO: clean

 ndof  = size(f,1);
 ntime = size(f,2);

 u = zeros( ndof, ntime );
 v = zeros( ndof, ntime );
 a = zeros( ndof, ntime );
 
 up = u0; vp = v0;
 ap = M\(f(:,1)-C*vp-K*up); 

 fatK = K + 1/(beta*dt^2)*M + gamma/(beta*dt)*C;
 
 invert = 0;
% if fatK1 != false
%    invert = 1;
 if 2*ntime > ndof % It's preferable to invert fatK
    fatK1 = inv(fatK);
    invert = 1;
 else % Pre-factorize
    [Rr,p,Qq] = chol (fatK, "vector");
    % TODO: control on p (should be zero for fatK to be sym. def. pos.)
 end
 
 doIneq = false;
 if isfield(ineq,"C") % Recover parameters for Uzawa
    Ci = ineq.C;  % C u <= d
    di = ineq.d;
    ki = ineq.k;
    ni = ineq.n;
    
    nineq  = size(Ci,1); % nb of inequalities to ensure
    doIneq = true; % Just a flag
 else
    ni = 1;
    nineq = 1;
 end
 
 u(:,1) = u0; v(:,1) = v0; a(:,1) = ap;
 
 for i=2:ntime
    b = ( f(:,i) + ...
         C*( gamma/(beta*dt)*up + (gamma/beta-1)*vp + dt/2*(gamma/beta-2)*ap ) + ...
         M*( 1/(beta*dt^2)*up + 1/(beta*dt)*vp + (1/(2*beta)-1)*ap ) );
 
    floc = zeros(nineq,1); % Local Uzawa force
    fu = zeros(ndof,1); % Additionnal Uzawa force (will stay zero in case there is no inequality)
    for iter = 1:ni % Iterations for Uzawa : only one will be done in case there is no inequality
       if invert
          u(:,i) = fatK1*(b+fu);
       else
          %u(Qq,i) = Uu \ ( Ll \ (b(Pp)) );
          u(Qq,i) = Rr \ ( Rr' \ (b(Qq)+fu(Qq)) ); % TODO: is it optimal ?
       end
       
       if doIneq % Update Uzawa RHS. TODO: stopping criterion
          res = Ci*u(:,i) - di; % This should be <=0
          floc = floc - ki*res; % Retroactive force
          floc = .5*(floc-abs(floc)); % floc is necessarly negative
          fu = Ci'*floc; % Uzawa additionnal force
          
          %if iter==1, norm(res(find(res>0))), end
       end
    end
    %norm(res(find(res>0))) % TODO: handle residual. At least a warning
               
    a(:,i) = 1/(beta*dt^2) * ( u(:,i) - up - dt*vp - dt^2/2*(1-2*beta)*ap );
    v(:,i) = vp + dt*( (1-gamma)*ap + gamma*a(:,i) );
    
    up = u(:,i); vp = v(:,i); ap = a(:,i);
 end

end
