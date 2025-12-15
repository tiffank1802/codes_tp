function [u,v,a] = newmark1stepMRHS (M, C, K, f, u0, v0, a0, dt, beta, gamma)
 % Computes one step of the Newmark argorithm
 nrhs = size(u0,2);

 fatK = K + 1/(beta*dt^2)*M + gamma/(beta*dt)*C;

 %b = ( f*ones(1,nrhs) + ...
 b = ( f + ...
       C*( gamma/(beta*dt)*u0 + (gamma/beta-1)*v0 + dt/2*(gamma/beta-1)*a0 ) + ...
       M*( 1/(beta*dt^2)*u0 + 1/(beta*dt)*v0 + (1/(2*beta)-1)*a0 ) );
       
 u = fatK\b; % Rem : it would be more effective to pre-choleskify this one outside the time loop, but less understandable
 a = 1/(2*beta*dt^2/2) * ( u - u0 - dt*v0 - dt^2/2*(1-2*beta)*a0 );
 v = v0 + dt*( (1-gamma)*a0 + gamma*a );
end
