%% Extended Kalman Filter (EKF) for 2D Car Navigation Example
% This script simulates a car moving in 2D and uses an EKF to estimate its state.
clear; clc; close all;

%% Simulation Parameters
dt = 0.1;               % Time step [s]
T = 20;                 % Total simulation time [s]
t = 0:dt:T;           % Time vector

% Vehicle parameters
L = 2.5;                % Wheelbase length [m]

% Process noise covariance (assumed noise on state propagation)
Q = diag([0.1, 0.1, 0.05, 0.1].^2);

% Measurement noise covariance (for GPS: measuring x and y)
R = diag([0.5, 0.5].^2);
% R = diag([0,0]);

%% Initial Conditions
% True state: [x; y; theta; v]
x_true = [0; 0; 0; 5];
% Initial state estimate
x_est = [0; 0; 0; 5];
% Initial estimation error covariance
P = diag([1, 1, 0.1, 1].^2);

% Preallocate arrays for storing history for plotting
x_true_history = zeros(4, length(t));
x_est_history = zeros(4, length(t));

%% Control Inputs (constant for simplicity)
a = 0.2;         % Acceleration [m/s^2]
delta = 0.05;    % Steering angle [radians]

%% EKF Simulation Loop
for k = 1:length(t)
    %% True State Propagation
    theta = x_true(3);
    v = x_true(4);
    % Update true state using the non-linear kinematic bicycle model
    x_true(1) = x_true(1) + v * cos(theta) * dt;
    x_true(2) = x_true(2) + v * sin(theta) * dt;
    x_true(3) = x_true(3) + (v / L) * tan(delta) * dt;
    x_true(4) = x_true(4) + a * dt;
    
    % Simulated GPS measurement (position only) with added noise
    z = [x_true(1); x_true(2)] + sqrt(R) * randn(2, 1);
%     z = [x_true(1); x_true(2)];
    
    %% EKF Prediction Step
    % Predict the next state using the same non-linear model
    theta_est = x_est(3);
    v_est = x_est(4);
    x_pred = zeros(4, 1);
    x_pred(1) = x_est(1) + v_est * cos(theta_est) * dt;
    x_pred(2) = x_est(2) + v_est * sin(theta_est) * dt;
    x_pred(3) = x_est(3) + (v_est / L) * tan(delta) * dt;
    x_pred(4) = x_est(4) + a * dt;
    
    % Compute the Jacobian F of the process model with respect to the state
    F = eye(4);
    F(1,3) = -v_est * sin(theta_est) * dt;
    F(1,4) = cos(theta_est) * dt;
    F(2,3) = v_est * cos(theta_est) * dt;
    F(2,4) = sin(theta_est) * dt;
    F(3,4) = (1 / L) * tan(delta) * dt;
    
    % Predict the error covariance
    P_pred = F * P * F' + Q;
    
    %% EKF Update Step
    % Measurement model: z = h(x) = [x; y]
    z_pred = [x_pred(1); x_pred(2)];
    % Measurement Jacobian H
    H = [1 0 0 0;
         0 1 0 0];
     
    % Innovation covariance
    S = H * P_pred * H' + R;
    % Kalman gain
    K = P_pred * H' / S;
    
    % Update the state estimate with the measurement
    x_est = x_pred + K * (z - z_pred);
    % Update the error covariance
    P = (eye(4) - K * H) * P_pred;
    
    % Save states for plotting
    x_true_history(:, k) = x_true;
    x_est_history(:, k) = x_est;
end

%% Plotting the Results
figure;
plot(x_true_history(1, :), x_true_history(2, :), 'b-', 'LineWidth', 2); hold on;
plot(x_est_history(1, :), x_est_history(2, :), 'r--', 'LineWidth', 2);
xlabel('X Position [m]');
ylabel('Y Position [m]');
legend('True Path', 'Estimated Path');
title('2D Car Navigation using EKF');

figure;
subplot(2,2,1);
plot(t, x_true_history(1, :), 'b-', t, x_est_history(1, :), 'r--');
xlabel('Time [s]'); ylabel('X Position [m]');
title('X Position');

subplot(2,2,2);
plot(t, x_true_history(2, :), 'b-', t, x_est_history(2, :), 'r--');
xlabel('Time [s]'); ylabel('Y Position [m]');
title('Y Position');

subplot(2,2,3);
plot(t, x_true_history(3, :), 'b-', t, x_est_history(3, :), 'r--');
xlabel('Time [s]'); ylabel('Heading \theta [rad]');
title('Heading');

subplot(2,2,4);
plot(t, x_true_history(4, :), 'b-', t, x_est_history(4, :), 'r--');
xlabel('Time [s]'); ylabel('Velocity [m/s]');
title('Velocity');

sgtitle('EKF State Estimation vs True State');
