%% 2D Localization with Sensor Noise and Kalman Filter
clc; clear;close all
% Simulation parameters
dt = 0.1;               % time step (seconds)
T = 20;                 % total simulation time (seconds)
t = 0:dt:T;
nSteps = length(t);

% True initial state: [x; y; vx; vy]
x_true = zeros(4, nSteps);
x_true(:,1) = [0; 0; 1; 0.5];

% Constant velocity motion model (state transition)
A = [1 0 dt 0;
     0 1 0 dt;
     0 0 1  0;
     0 0 0  1];

% Process noise (simulate small acceleration effects)
sigma_acc = 0.2;
Q = [dt^4/4,    0, dt^3/2,    0;
         0, dt^4/4,    0, dt^3/2;
     dt^3/2,    0,   dt^2,    0;
         0, dt^3/2,    0,   dt^2] * sigma_acc^2;

% Simulate true trajectory with process noise
for k = 2:nSteps
    % Adding process noise to the state propagation
    x_true(:,k) = A * x_true(:,k-1) + mvnrnd([0;0;0;0], Q)';
end

% Measurement model: we only measure the position [x; y]
H = [1 0 0 0;
     0 1 0 0];
sigma_meas = 0.5;    % standard deviation of the measurement noise
R = sigma_meas^2 * eye(2);

% Simulate noisy measurements
z = zeros(2, nSteps);
for k = 1:nSteps
    z(:,k) = H * x_true(:,k) + sigma_meas * randn(2,1);
end

% Kalman Filter Initialization
x_est = zeros(4, nSteps);      % state estimates
P = eye(4);                    % initial covariance estimate
x_est(:,1) = [0; 0; 1; 0.5];     % initial state estimate

% Kalman Filter Loop
for k = 2:nSteps
    % Prediction step
    x_pred = A * x_est(:,k-1);
    P_pred = A * P * A' + Q;
    
    % Measurement update
    K = P_pred * H' / (H * P_pred * H' + R);  % Kalman gain
    x_est(:,k) = x_pred + K * (z(:,k) - H * x_pred);
    P = (eye(4) - K * H) * P_pred;
    plot(x_est(1,:), x_est(2,:), 'g--', 'LineWidth', 2);
end

% Plotting the results
figure;
plot(x_true(1,:), x_true(2,:), 'b-', 'LineWidth', 2); hold on;
plot(z(1,:), z(2,:), 'rx', 'MarkerSize', 4);
plot(x_est(1,:), x_est(2,:), 'g--', 'LineWidth', 2);
legend('True Trajectory', 'Noisy Measurements', 'Estimated Trajectory');
xlabel('X Position');
ylabel('Y Position');
title('2D Localization with Sensor Noise and Kalman Filter');
grid on;
