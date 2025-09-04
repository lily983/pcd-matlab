% Plant Model (example)
Ts = 0.1; % Sample time
A = [1 0.1; 0 1];
B = [0.005; 0.005];
C = [1 0; 0 1];
D = [0; 0];
plant = ss(A, B, C, D, Ts);

% MPC Controller Design
mpcobj = mpc(plant);
mpcobj.PredictionHorizon = 10;
mpcobj.ControlHorizon = 3;
mpcobj.Weights.Output = [1, 1]; % Weight on x and y errors
mpcobj.Weights.Input = [0.1, 0.1]; % Weight on u1 and u2 (control inputs)

% Simulate
u = [1; 1]; % Initial control input
y = [0; 0]; % Initial state
r = [1; 1]; % Reference trajectory
[y, u] = mpcsim(mpcobj, u, y, r);

% Plot Results
plot(mpcobj, y, r, u);