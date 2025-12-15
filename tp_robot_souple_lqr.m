% 22/11/2022 : TP pilotage d'un robot souple via Kalman et le contrôle optimal
% Exercise 4.2.2: LQR Control Implementation

close all;
clear all;

nstep = 1000; % Nb time steps
delta = .01; % Time step
q     = .0001; % Error on model
r     = .0001; % Error on measurement
T      = 3; % command Period

mu    = 1e-1; % Optimal control regularity parameter

beta = .25;
gamma = .5; % Newmark parameters

L     = 500; % Beam length
a     = 5;
b     = 5; % Beam section dimensions
E     = 70000; % Young modulus
rho   = 2700e-9; % Mass
cy    = 1e-4; % Damping coefficient
nx    = 10; % total nb of beam elements
nA    = floor(nx/2); % nb of beam elements to point A

S  = a*b;   % Section
Ix = a*b^3/12; % Moment

dx = L/nx;
time = delta:delta:delta*nstep;
time0 = [0,time];
ixe  = dx:dx:L;

isMatlab = true;
vseed1 = 0; vseed2 = 1; vseed3 = 2;

% Get indices of points
npt = nx + 1; ndof = 2*npt; % nb dofs before BCs
ndo1 = ndof-2; % nb of dofs after BCs
induA = 2*nA-1;
induB = 2*nx-1;

% Perturbations
fpert = generateNoiseTemporal( time0, 1, q, vseed1, isMatlab ); % Perturbation force
mpert = generateNoiseTemporal( time0, 1, r, vseed2, isMatlab ); % Measure perturbation

% Operators and BCs
fA_vec = zeros(1,size(time0,2));
ffull = zeros(ndof,nstep+1);
[Kfull,Cfull,Mfull] = neb_beam_matrices( nx, dx, E, Ix, rho, S, cy ); % Stiffness & Co.
K = Kfull; K(:,[1,2]) = []; K([1,2],:) = [];
C = Cfull; C(:,[1,2]) = []; C([1,2],:) = [];
M = Mfull; M(:,[1,2]) = []; M([1,2],:) = [];
f = ffull; f([1,2],:) = [];

% Initial conditions
u0 = zeros(ndo1,1); v0 = zeros(ndo1,1);

%% LQR Design
% State space formulation: x = [u; v]
% M*a + C*v + K*u = F
% a = -M\K*u - M\C*v + M\F
% x_dot = [0 I; -M\K -M\C] * x + [0; M\Bf] * fA

n_state = 2 * ndo1;
Minv = full(M) \ eye(size(M)); % Precompute inverse (or solve linear system)
MinvK = M \ K;
MinvC = M \ C;

A = [zeros(ndo1), eye(ndo1);
     -MinvK,     -MinvC];

% Input matrix B
B_force = zeros(ndo1, 1);
B_force(induA) = 1; % Force applied at induA
B = [zeros(ndo1, 1);
     M \ B_force];

% LQR Weights
% Q penalizes state deviation. Let's penalize the tip position error.
% Output y = C_out * x. We want y -> ref.
% C_out selects induB from u.
C_out = zeros(1, n_state);
C_out(induB) = 1;

% Q = C_out' * C_out * 100; % Penalize tip position
% Alternatively, penalize all states to keep them small (regulator)
Q = eye(n_state); 
% Or maybe just the positions?
% Q = diag([ones(ndo1,1); zeros(ndo1,1)]);

% Let's use a simple identity for now, or the one implied by the problem.
% Given 'mu' is the only parameter, maybe Q is fixed (e.g. Identity or C'C).
% Let's try Q = C_out'*C_out + epsilon*I to ensure observability?
% Let's stick to Q = eye(n_state) for generic stabilization, 
% or better, penalize the error we want to minimize.
% The PID minimized (1 - u_tip).
% So let's put high weight on u_tip.
Q = zeros(n_state);
Q(induB, induB) = 1e4; % High weight on tip position
% Add small regularization on everything else to avoid singular Q if needed (though lqr handles it)
Q = Q + 1e-3 * eye(n_state);

R = mu;

% Compute LQR gain
fprintf('Computing LQR gain...\n');
[K_lqr, S_lqr, E_lqr] = lqr(A, B, Q, R);

% Precompute Feedforward for tracking
% We want y = C_out * x = 1 (Step reference)
% At steady state:
% [A B; C_out 0] * [x_ss; u_ss] = [0; 1]
% Solve for x_ss, u_ss
BigMat = [A, B; C_out, 0];
RHS = [zeros(n_state, 1); 1];
sol = BigMat \ RHS;
x_ss = sol(1:n_state);
u_ss = sol(n_state+1);

fprintf('Steady state control u_ss: %f\n', u_ss);

%% Simulation
u = zeros(ndo1,nstep+1);
v = zeros(ndo1,nstep+1);
a = zeros(ndo1,nstep+1);
a0 = M \ (f(:,1)-C*v0-K*u0);
a(:,1) = a0;

% Initialize state estimate (perfect knowledge for now)
x_curr = [u0; v0];

for step = 2:nstep+1
    % Current state (perfectly observed)
    x_curr = [u(:,step-1); v(:,step-1)];
    
    % LQR Control Law with Feedforward
    % u_ctrl = -K * (x - x_ss) + u_ss
    f_ctrl = -K_lqr * (x_curr - x_ss) + u_ss;
    
    % Add perturbation
    fA_val = f_ctrl ;%+ fpert(step);
    
    % Store force
    f(induA,step) = fA_val;
    fA_vec(step) = fA_val;
    
    % Newmark integration
    [u(:,step), v(:,step), a(:,step)] = newmark1stepMRHS(...
        M, C, K, f(:,step), u(:,step-1), v(:,step-1), a(:,step-1), delta, beta, gamma);
end

%% Plotting
figure; hold on;
plot(ixe, u(1:2:end-1,end),'Color','blue','linewidth',3);
plot(ixe, u(2:2:end,end)  ,'Color','red','linewidth',3);
h = legend('u','\theta'); title('U_{end} (LQR)');
xlabel('x'); ylabel('v/\theta');
set(gca,'FontSize',16);
set(h,'FontSize',16);

figure; hold on;
plot(time0, u(induB,:),'Color','blue','linewidth',3);
plot(time0, fA_vec,'Color','red','linewidth',3);
yline(1, '--k', 'Reference');
h = legend('u_{tip}','fA', 'Ref'); title('U_{time} (LQR)');
xlabel('t'); ylabel('v');
set(gca,'FontSize',16);
set(h,'FontSize',16);
