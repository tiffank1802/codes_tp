function fpert = generateNoiseTemporal( time, tlength, q, vseed, isMatlab )
   % Generates correlated noise
   
   nstep = max(size(time));
   
   if isMatlab, rng(vseed); else randn('state',vseed); end % Set random seed
   seed  = randn(nstep,1);                  % basic noise
   
   Corr = exp(-((time-time')).^2/tlength^2);  % (square root of) Correlation matrix
   fpert = Corr*seed;                       % correlated noise
   
   famp = sum(fpert.^2)/nstep;
   fpert = fpert * sqrt(q/famp); % Tune to get the right amplitude
end
