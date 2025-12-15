function [fop, ress] = forthodir( m, c, k, co, u0, v0, niter, nopt, stagEps, delta, beta, gamma, mu, picture, indices, action ) % TODO: more arguments?
   % Optimize loading to follow a given prescripted displacement Using an ORTHODIR algorithm
   
   % indices : components of u that must respect the consign
   % action : components of f that can be controlled

   ndof = size(u0,1);
   if ~exist("indices")
      indices = 1:ndof; % Use all components
   end
   if ~exist("action")
      action = 1:ndof; % Use all components
   end
   nact = max(size(action)); % Nb of actuators

    fop   = zeros(ndof,niter); % TODO: allow non-zero (mainly for non-causal controller)
   res   = zeros(nopt,1);
   Delta = zeros(nopt,1);
    dirs  = zeros(nact*niter,nopt); % Research direction storage
    Ads   = zeros(nact*niter,nopt);
   
    diff   = zeros(ndof, niter); % Will store differences
    dirdir = zeros(ndof, niter);

   % Determine RHS
   ai(:,1) = - m \ (c*v0 + k*u0); % Initial acceleration (should actually probably be inside Newmark instead)
    [ui,vi,ai]  = Newmark (m, c, k, zeros(ndof,niter), u0, v0, ai(:,1), delta, beta, gamma); % Inverse dynamics % TODO: re-code with Newmark2 instead
    diff(indices,:) = co(indices,:) - ui(indices,:);
   laa(:,1) = m \ diff(:,end); % Final value
   [lau,lav,laa] = Newmark (m, c, k, diff(:,end:-1:1), 0, 0, laa(:,1), delta, beta, gamma); % Inverse dynamics
   laulau = lau(action,end:-1:1); % Right hand side
   rhs = laulau(:); % Reshape as a vector

   % Residual b-Ax and search direction
   resv = rhs; % Initialization = 0
   res0 = norm(resv);
   
   % Test to avoid singularity
   if res0 < stagEps
      disp('zero seems to be the optimal solution');
      return;
   end

   % Loop
   if picture, h = waitbar( 0, 'Optimization' ); end
   for iter=1:nopt
      % New directions
      dir = resv;
      
      % Determine A*dir
       dirdir(action,:) = reshape(dir,[nact,niter]); % Reshape as a multi-vector
      ai(:,1) = m \ dirdir(:,1); % Initial acceleration
      [ui,vi,ai]  = Newmark (m, c, k, dirdir, 0, 0, ai(:,1), delta, beta, gamma); % Inverse dynamics
      diff(indices,:) = (-ui(indices,:));
      laa(:,1) = m \ diff(:,end); % Final value
      [lau,lav,laa] = Newmark (m, c, k, diff(:,end:-1:1), 0, 0, laa(:,1), delta, beta, gamma); % Inverse dynamics
      laulau = lau(action,end:-1:1);
      Ad = mu*dir - laulau(:);
      
      % Orthogonalize
      for j=1:iter-1
         phiij = Ads(:,j)'*Ad; betaij = phiij/Delta(j);
         dir = dir - betaij*dirs(:,j);
         Ad  = Ad  - betaij*Ads(:,j);
      end
      
      % Store directions
      dirs(:,iter) = dir;
      Ads(:,iter)  = Ad;

      % Step
      Delta(iter) = Ad'*Ad; gammai = Ad'*resv; alphai = gammai/Delta(iter);
      
      % Apply
       dirdir(action,:) = reshape(dir,[nact,niter]);
      fop  = fop + alphai*dirdir;
      resv = resv - alphai*Ad;
      res(iter) = norm(resv);
      
      % Convergence test
      if res(iter)/res0 < stagEps
         resi = res(iter);
         break;
      end

      if picture, waitbar(iter/nopt,h); end
   end
   if picture, delete(h); end

   ress = res(1:iter);
end
