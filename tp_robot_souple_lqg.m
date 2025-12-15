% 22/11/2022 : TP pilotage d'un robot souple via Kalman et le contrôle optimal
% Part 4: LQG Control (LQR + Kalman Filter)

close all;
clear all;

% Parameters
nstep = 1000; % Nb time steps
delta = .01; % Time step
q     = .0001; % Error on model (Process noise variance)
r     = .0001; % Error on measurement (Measurement noise variance)
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
[Kfull,Cfull,Mfull] = neb_beam_matrices( nx, dx, E, Ix, rho, S, cy ); % Stiffness & Co.
K = Kfull; K(:,[1,2]) = []; K([1,2],:) = [];
C = Cfull; C(:,[1,2]) = []; C([1,2],:) = [];
M = Mfull; M(:,[1,2]) = []; M([1,2],:) = [];

% Initial conditions
u0 = zeros(ndo1,1); v0 = zeros(ndo1,1);

%% State Space Model
% x = [u; v]
n_state = 2 * ndo1;
Minv = full(M) \ eye(size(M));
MinvK = M \ K;
MinvC = M \ C;

A_c = [zeros(ndo1), eye(ndo1);
       -MinvK,     -MinvC];

% Input matrix B (Force at A)
B_force = zeros(ndo1, 1);
B_force(induA) = 1; 
B_c = [zeros(ndo1, 1);
       M \ B_force];

% Output matrix C (Tip position)
C_out = zeros(1, n_state);
C_out(induB) = 1;

% Discretize System
sys_c = ss(A_c, B_c, C_out, 0);
sys_d = c2d(sys_c, delta, 'zoh');
[Ad, Bd, Cd, Dd] = ssdata(sys_d);

%% 1. LQR Design (Discrete)
% Penalize tip position error
Q_lqr = zeros(n_state);
Q_lqr(induB, induB) = 1e4; 
Q_lqr = Q_lqr + 1e-3 * eye(n_state); % Regularization
R_lqr = mu;

fprintf('Computing Discrete LQR gain...\n');
[K_lqr, S_lqr, E_lqr] = dlqr(Ad, Bd, Q_lqr, R_lqr);

% Feedforward for tracking
% [Ad-I  Bd; Cd 0] * [x_ss; u_ss] = [0; 1]
BigMat = [Ad - eye(n_state), Bd; Cd, 0];
RHS = [zeros(n_state, 1); 1];
sol = BigMat \ RHS;
x_ss = sol(1:n_state);
u_ss = sol(n_state+1);

%% 2. Kalman Filter Design
% Model: x[k+1] = Ad x[k] + Bd (u[k] + w[k])
%        y[k]   = Cd x[k] + v[k]
% w ~ N(0, q), v ~ N(0, r)
% Process noise covariance matrix Q_kal
% Since w is scalar scaling B_d:
Q_kal = q; 
R_kal = r;

% Compute Kalman Gain
% We use dlqe. G = Bd.
fprintf('Computing Kalman gain...\n');
[L_kal, P_kal, Z_kal] = dlqe(Ad, Bd, Cd, Q_kal, R_kal);

%% Simulation Loop
u = zeros(ndo1,nstep+1);
v = zeros(ndo1,nstep+1);
a = zeros(ndo1,nstep+1);
f = zeros(ndo1,nstep+1);

% Initial state estimate
x_est = zeros(n_state, 1); 
x_est_hist = zeros(n_state, nstep+1);

% Initial acceleration for Newmark
a0 = M \ (-C*v0 - K*u0); 
a(:,1) = a0;

% Storage for plotting
y_meas_hist = zeros(1, nstep+1);
u_ctrl_hist = zeros(1, nstep+1);

for step = 2:nstep+1
    % 1. Measurement (from previous step real state)
    % y[k-1] = C * x[k-1] + v[k-1]
    y_meas = u(induB, step-1) + mpert(step-1);
    y_meas_hist(step-1) = y_meas;
    
    % 2. Kalman Prediction & Update
    % We are at step 'step'. We want to compute control for 'step'.
    % We have measurement up to 'step-1'.
    % Standard digital control loop:
    % Measure y[k]. Update estimate x[k|k]. Compute u[k]. Apply u[k].
    % Here 'step' corresponds to time index k.
    
    % Prediction step (Time Update) from k-1 to k
    if step == 2
        x_pred = x_est; % Initial guess
    else
        % x_pred = Ad * x_est_prev + Bd * u_prev
        x_pred = Ad * x_est + Bd * u_ctrl_hist(step-1);
    end
    
    % Measurement Update (Correction)
    % Innovation: y[k] - Cd * x_pred
    % Note: In this loop structure, we are measuring the state at 'step-1' (which is the result of previous integration).
    % So let's align indices carefully.
    % u(:, step-1) is state at t = (step-2)*delta.
    % Wait, time0 = [0, delta, ...].
    % step=1 -> t=0. step=2 -> t=delta.
    % u(:,1) is initial condition at t=0.
    % u(:,2) is computed from u(:,1) via Newmark.
    
    % At start of loop step=2 (t=delta):
    % We have u(:,1) (t=0).
    % We want to compute f(:,2) to go from 1 to 2?
    % Newmark computes u(:,step) using f(:,step).
    % So f(:,step) is the force applied at t=delta.
    % We need to decide f(:,step) based on info available.
    % Info available: Measurement at t=0? Or t=delta?
    % Usually, we measure y(t), compute u(t), apply u(t).
    % So at step=2, we measure y(step-1) (t=0)?
    % If we measure y(t=0), we estimate x(0), compute u(0)?
    % But u(:,1) is already fixed (IC).
    % We are computing f(:,step) which drives the system to step.
    % Actually Newmark: M a_n+1 + C v_n+1 + K u_n+1 = F_n+1
    % F_n+1 is force at t_{n+1}.
    % So at step=2, we need F at t=delta.
    % We can measure y at t=delta? No, y at t=delta depends on u at t=delta which depends on F.
    % So we must use y at t=0 (step-1).
    % So we use estimate x_est(step-1) (at t=0) projected to t=delta?
    % Or we use x_est(step-1) to compute u(step-1)?
    % But F is applied at step.
    
    % Let's assume:
    % Measure y[k-1].
    % Update estimate x[k-1|k-1].
    % Predict x[k|k-1] = Ad * x[k-1|k-1] + Bd * u[k-1].
    % Compute u[k] = -K * x[k|k-1].
    
    % Measurement at step-1
    y_k_minus_1 = u(induB, step-1) + mpert(step-1);
    
    % If step=2 (k=1), we measure y[0].
    % x_est initialized to 0 (at k=0).
    % Correction at k-1:
    if step == 2
        x_est_prev = x_est; % x[0|0] approx
        % Or better:
        % x_est_prev = x_est + L_kal * (y_k_minus_1 - Cd * x_est);
    else
        x_est_prev = x_est; % This is x[k-1|k-1] from previous loop
    end
    
    % Correction (Update) for step-1
    % Actually, let's keep track of x_pred and x_corr.
    % Let x_est variable hold x[k-1|k-1].
    if step == 2
        x_pred_k_minus_1 = zeros(n_state,1); % x[0|-1]
        u_prev = 0;
    else
        % From previous iteration, we have x_est (which was x[k-2|k-2])
        % and u_ctrl_hist(step-2).
        % Wait, let's simplify.
        
        % Let's maintain x_est as x[k-1|k-1].
        % At step=2 (target t=delta), we have measurement at step-1 (t=0).
        % We update x[0|0] using y[0].
        % Then predict x[1|0] using u[0].
        % Then compute u[1] using x[1|0].
    end
    
    % Correct logic:
    % 1. Measure y(t_{k-1}).
    % 2. Update estimate: x_{k-1|k-1} = x_{k-1|k-2} + L(y_{k-1} - C x_{k-1|k-2})
    % 3. Predict next: x_{k|k-1} = A x_{k-1|k-1} + B u_{k-1}
    % 4. Control: u_k = -K x_{k|k-1} + u_ss
    
    % Need to store x_pred (x_{k|k-1}) between steps.
    if step == 2
        x_pred = zeros(n_state, 1); % x[0|-1] (Initial guess)
        u_prev = 0; % u[-1]
    else
        x_pred = x_next_pred; % Retrieved from previous step
        u_prev = u_ctrl_hist(step-1); % u[k-1] used in prediction? 
        % Wait, u_ctrl_hist(step-1) is u at t_{k-1}.
        % Yes.
    end
    
    % Measurement
    y_meas = u(induB, step-1) + mpert(step-1);
    
    % Update
    innovation = y_meas - Cd * x_pred;
    x_est = x_pred + L_kal * innovation; % x[k-1|k-1]
    
    % Save estimate
    x_est_hist(:, step-1) = x_est;
    
    % Control for CURRENT step (t_k)
    % We need x[k|k-1] to compute u[k].
    % But we haven't computed u[k-1] yet? 
    % No, u[k-1] was computed in previous step.
    % Wait, at step=2, we are computing f(:,2) (u[1]).
    % We used u[0] (f(:,1)) previously.
    % f(:,1) was 0 (or defined by IC).
    
    % Predict next state for control
    % x[k|k-1] = Ad * x[k-1|k-1] + Bd * u[k-1]
    % Here u[k-1] is the control applied at previous step.
    % For step=2 (k=1), u[k-1] = u[0] = f(:,1) at induA.
    u_applied_prev = f(induA, step-1);
    
    x_next_pred = Ad * x_est + Bd * u_applied_prev;
    
    % Compute Control u[k]
    % u[k] = -K * (x[k|k-1] - x_ss) + u_ss
    f_ctrl = -K_lqr * (x_next_pred - x_ss) + u_ss;
    
    % Apply to system
    % Add perturbation to the plant input
    fA_val = f_ctrl + fpert(step);
    
    f(induA, step) = fA_val;
    u_ctrl_hist(step) = f_ctrl;
    
    % Simulate Plant (Newmark)
    [u(:,step), v(:,step), a(:,step)] = newmark1stepMRHS(...
        M, C, K, f(:,step), u(:,step-1), v(:,step-1), a(:,step-1), delta, beta, gamma);
        
end

%% Plotting
figure('Name', 'LQG Control - Tip Position');
subplot(2,1,1);
hold on;
plot(time0, u(induB,:), 'b', 'LineWidth', 2);
plot(time0, x_est_hist(induB, :), 'g--', 'LineWidth', 2);
yline(1, 'k--');
legend('Real Tip Pos', 'Estimated Tip Pos', 'Reference');
title('Tip Position Tracking (LQG)');
grid on;

subplot(2,1,2);
plot(time0, f(induA,:), 'r');
title('Control Input (Force at A)');
xlabel('Time (s)');
grid on;

figure('Name', 'Estimation Error');
plot(time0, u(induB,:) - x_est_hist(induB, :));
title('Estimation Error at Tip');
xlabel('Time (s)');
grid on;

% Compare with LQR (if desired, but this is standalone)
