clc
clear 
close all


t_sim = 30;
t_sample = 0.03;

%% Reference generation
t_step = t_sim/5; % time for the step input
t_stepdown = t_sim*4/5;

% Time vector for simulation
t_vec = 0:t_sample:t_sim;

x1ref_before = zeros(find(t_step == t_vec),1)+1;
step_size = 3;
x1ref_step = step_size*ones(length(t_vec)-length(x1ref_before),1);
x1ref = [x1ref_before;x1ref_step];
idx_step_down = find(t_stepdown == t_vec);
x1ref(idx_step_down:end) = 0;


x2ref_before = zeros(find(t_step == t_vec),1)+1;
step_size = 2;
x2ref_after = step_size*ones(length(t_vec)-length(x2ref_before),1);
x2ref = [x2ref_before;x2ref_after];
idx_step_down = find(t_stepdown == t_vec);
x2ref(idx_step_down:end) = 0;


%% Plant (replace by linmod) 
% Mass in kg
m1 = 1;
m2 = 1.2;

% Spring constants 
k1 = 20;
k2 = 15;

% Damping coefficient
c1 = 3;
c2 = 2;

A = [0 1 0 0;
    -(k1+k2)/m1 -(c1+c2)/m1 k2/m1 c2/m1;
     0 0 0 1;
     k2/m2 c2/m2 -k2/m2 -c2/m2];

B = [0 0;
     1/m1 0;
     0 0;
     0 1/m2];

C = [1 0 0 0;
     0 0 1 0];

D = zeros(2,2);

sys = ss(A,B,C,D);

%% MPC preparation using three penalty method 
% Discretize the state-space model for simulation
sys_d = c2d(sys, t_sample);
Ad = sys_d.a;
Bd = sys_d.b;
Cd = sys_d.c;
Dd = sys_d.d;

[m, p] = size(Bd); % m states, p inputs
[n, ~] = size(Cd); % n outputs

% horizon sizes
Np = 10;
Nc = Np;

% Offline computation

% prediction matrix F G
[F,G] = predmat(Ad,Bd,Cd,Np); 

Qy  = diag([10 10]);     % penalize y, match n dimension
Qu  = diag([0 0]);   % penalize input magnitude, match p dimension
Qur = diag([1e-4 1e-4]);       % penalize input increments, match p dimension

% Cost weight
Qs   = kron(eye(Np), Qy);
Qus  = kron(eye(Np), Qu);
Qurs = kron(eye(Np), Qur);

I = eye(p);
M = zeros(p*Np);
c = zeros(p*Np,p);

% first block
M(1:p,1:p) = I;
c(1:p,:) = -I;

% remaining blocks
for i = 2:Np
    
    row = (i-1)*p + (1:p);
    
    M(row,(i-2)*p + (1:p)) = -I;
    M(row,(i-1)*p + (1:p)) = I;
    
end

%% offline QP prep

H = 2 * (G'*Qs*G + Qus + M'*Qurs*M);
fx = 2*G'*Qs*F;
fr = -2*G'*Qs;
fu = -2*M'*Qurs*c;


% Constrained MPC
lb = -500*ones(p*Np,1);
ub =  500*ones(p*Np,1);


%% KF offline preparation

% covariance matrix
Q = 0.01*eye(m); % model noise (tunning knob)
R = 0.0001*eye(n); % sensor noise (fixed) 

if rank(obsv(Ad,Cd)) ~= m
    disp('System not obsv able')
elseif rank(ctrb(A,B)) ~= m
    disp("System not ctrl able")
end

% sim("MPC_example_2mass_2024a.slx")