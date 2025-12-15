close all;
clear all;

%% Paramètres temporels
nstep = 1000; % Nb time steps
delta = .01; % Time step
q     = .0001; % Error on model (perturbation variance)
r     = .0001; % Error on measurement (non utilisé ici car pas de bruit)
T     = 3; % command Period

%% Paramètres PID
Kp    = 5e-1;
Ki    = 1e-1;
Kd    = 1e-2;

beta = .25;
gamma = .5; % Newmark parameters

%% Paramètres poutre
L     = 500; % Beam length
a     = 5;
b     = 5; % Beam section dimensions
E     = 70000; % Young modulus
rho   = 2700e-9; % Mass
cy    = 1e-4; % Damping coefficient
nx    = 10; % total nb of beam elements
% nA    = floor(nx/2); % point A (milieu poutre)
 nA    = floor(nx); % point A (à l'extremité de la poutre)


S  = a*b;   % Section
Ix = a*b^3/12; % Moment
dx = L/nx;
time = delta:delta:delta*nstep;
time0 = [0,time];
ixe  = dx:dx:L;

isMatlab = true;
vseed1 = 0; vseed2 = 1; vseed3 = 2;

% Indices
npt = nx + 1; ndof = 2*npt; % nb dofs before BCs
ndo1 = ndof-2; % nb of dofs after BCs
induA = 2*nA-1;
induB = 2*nx-1;

%% Perturbation externe
fpert = generateNoiseTemporal( time0, 1, q, vseed1, isMatlab ); % Perturbation force
mpert = generateNoiseTemporal( time0, 1, r, vseed2, isMatlab ); % Measure perturbation


%% Matrices EF
[Kfull,Cfull,Mfull] = neb_beam_matrices( nx, dx, E, Ix, rho, S, cy );
K = Kfull; K(:,[1,2]) = []; K([1,2],:) = [];
C = Cfull; C(:,[1,2]) = []; C([1,2],:) = [];
M = Mfull; M(:,[1,2]) = []; M([1,2],:) = [];

%% Conditions initiales
u0 = zeros(ndo1,1); 
v0 = zeros(ndo1,1);

%% --- PID control avec perturbation et sans bruit ---
u = zeros(ndo1,nstep+1);
v = zeros(ndo1,nstep+1);
a = zeros(ndo1,nstep+1);

% Initialisation
u(:,1) = u0;
v(:,1) = v0;
a(:,1) = M \ (zeros(ndo1,1) - C*v0 - K*u0);

% Consigne Heaviside
t0 = 1.0;
u_ref = double(time0 >= t0);

% PID memory
integ_e = 0;
prev_e  = 0;

for k = 2:nstep+1
    % --- Mesure sans bruit ---
    uB_meas = u(induB, k-1)+mpert(k);
    
    % --- Erreur ---
    e = u_ref(k) - uB_meas;
    integ_e = integ_e + e*delta;
    de = (e - prev_e)/delta;
    
    % --- Commande PID ---
    Fpid = Kp*e + Ki*integ_e + Kd*de;
    
    % --- Force appliquée ---
    f_k = zeros(ndo1,1);
    f_k(induB) = Fpid;              % commande PID au point B
    f_k(induA) = f_k(induA) + fpert(k); % perturbation externe au point A
    
    % --- Newmark step ---
    [u1,v1,a1] = newmark1stepMRHS(M, C, K, f_k, ...
        u(:,k-1), v(:,k-1), a(:,k-1), ...
        delta, beta, gamma);
    
    u(:,k) = u1; v(:,k) = v1; a(:,k) = a1;
    prev_e = e;
end

%% --- Tracés PID avec perturbation ---
figure; hold on;
plot(time0, u(induB,:), 'b', 'LineWidth', 2);
plot(time0, u_ref, 'r--', 'LineWidth', 2);
legend('u_B(t)','Consigne');
xlabel('t (s)'); ylabel('Déplacement u_B');
title('Boucle fermée PID avec perturbation (sans bruit de mesure)');
set(gca,'FontSize',14);

figure; hold on;
plot(ixe, u(1:2:end-1,end),'Color','blue','linewidth',3);
plot(ixe, u(2:2:end,end) ,'Color','red','linewidth',3);
legend('u','\theta');
title('Déformée finale (PID avec perturbation seule)');
xlabel('x'); ylabel('u / \theta');
set(gca,'FontSize',14);
