% 22/11/2022 : TP pilotage d'un robot souple via Kalman et le contrôle optimal

close all;
clear all;

nstep = 1000; % Nb time steps
delta = .01; % Time step
q     = .0001; % Error on model
r     = .0001; % Error on measurement
T     = 3; % command Period

Kp    = 5e-1; % PID parameters
Ki    = 1e-1;%1e-2;%5e-2;
Kd    = 1e-2;%2e-2;%1e-1;%1e-2;%1e-5;
mu    = 1e-1; % Optimal control regularity parameter

beta = .25;
gamma = .5; % Newmark parameters

L     = 500; % Beam length
a     = 5;
b     = 5; % Beam section dimensions
E     = 70000; % Young modulus
rho   = 2700e-9; % Mass
cy    = 1e-4; % Damping coefficient % TODO: see if we keep that
nx    = 10; % total nb of beam elements
nA    = floor(nx/2);%floor(9*nx/10); % nb of beam elements to point A (watch out nx has to be multiple of 2)

S  = a*b;   % Section
Ix = a*b^3/12; % Moment

dx = L/nx;
time = delta:delta:delta*nstep;
time0 = [0,time];
ixe  = dx:dx:L;

isMatlab = true;
vseed1 = 0; vseed2 = 1; vseed3 = 2;

%%% Uncomment if running on OCTAVE (and comment if running on MATLAB).
%isMatlab = false;
%vseed1 = load('randomstate_head/randomState.mat');  vseed1 = vseed1.v;
%vseed2 = load('randomstate_head/randomState2.mat'); vseed2 = vseed2.v;
%vseed3 = load('randomstate_head/randomState3.mat'); vseed3 = vseed3.v;
%%% End Uncomment

% Get indices of points
npt = nx + 1; ndof = 2*npt; % nb dofs before BCs
ndo1 = ndof-2; % nb of dofs after BCs
induA = 2*nA-1;
induB = 2*nx-1;

% Perturbations
fpert = generateNoiseTemporal( time0, 1, q, vseed1, isMatlab ); % Perturbation force
mpert = generateNoiseTemporal( time0, 1, r, vseed2, isMatlab ); % Measure perturbation

% Operators and BCs
fA = sin(2*pi*time0/T);
ffull = zeros(ndof,nstep+1);
[Kfull,Cfull,Mfull] = neb_beam_matrices( nx, dx, E, Ix, rho, S, cy ); % Stiffness & Co.
K = Kfull; K(:,[1,2]) = []; K([1,2],:) = [];
C = Cfull; C(:,[1,2]) = []; C([1,2],:) = [];
M = Mfull; M(:,[1,2]) = []; M([1,2],:) = [];
f = ffull; f([1,2],:) = [];

% Set prescribed displacement at point A
co = zeros(ndo1, nstep+1);
co(induA,:) = fA; % desired displacement signal at point A
indices = induA;
action = induA;
niter = nstep + 1;
nopt = 100;
stagEps = 1e-6;
picture = 0;
addpath('old');

% Initial conditions
u0 = zeros(ndo1,1); v0 = zeros(ndo1,1);

% Direct resolution (shortcut)
% [u,v,a] = Newmark2N ( M, C, K, f, u0, v0, delta, beta, gamma );
% [f, ~] = forthodir(M, C, K, co, u0, v0, niter, nopt, stagEps, delta, beta, gamma, mu, picture, indices, action);

%% Direct resolution (step-by-step way)
u = zeros(ndo1,nstep+1);
v = zeros(ndo1,nstep+1);
a = zeros(ndo1,nstep+1);
a0 = M \ (f(:,1)-C*v0-K*u0); % initial acceleration
a(:,1) = a0;
% Consigne Heaviside
t0 = 1.0;
u_ref = double(time0 >= t0);

% PID memory
integ_e = 0;
prev_e  = 0;
for step = 2:nstep+1
  % --- Mesure sans bruit ---
    uB_meas = u(induB, step-1)+mpert(step);
    
    % --- Erreur ---
    e = u_ref(step) - uB_meas;
    integ_e = integ_e + e*delta;
    de = (e - prev_e)/delta;
    
    % --- Commande PID ---
    Fpid = Kp*e + Ki*integ_e + Kd*de;
    
    % --- Force appliquée ---
    f_k = zeros(ndo1,1);
    f_k(induB) = Fpid;              % commande PID au point B
    f_k(induA) = f_k(induA) + fpert(step); % perturbation externe au point A
    
    % --- Newmark step ---
    [u1,v1,a1] = newmark1stepMRHS(M, C, K, f_k, ...
        u(:,step-1), v(:,step-1), a(:,step-1), ...
        delta, beta, gamma);
    
    u(:,step) = u1; v(:,step) = v1; a(:,step) = a1;
    prev_e = e;
    end

% Plot some stuff
figure; hold on;
plot(ixe, u(1:2:end-1,end),'Color','blue','linewidth',3);
plot(ixe, u(2:2:end,end)  ,'Color','red','linewidth',3);
h = legend('u','\theta'); title('U_{end}');
xlabel('x'); ylabel('v/\theta');
set(gca,'FontSize',16);
set(h,'FontSize',16);

figure; hold on;
plot(time0, u(induB,:),'Color','blue','linewidth',3);
plot(time0, fA,'Color','red','linewidth',3);
h = legend('u','fA'); title('U_{time}');
xlabel('t'); ylabel('v');
set(gca,'FontSize',16);
set(h,'FontSize',16);

