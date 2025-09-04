%% EKF for Robot and Obstacle with Process Noise and Extraction of Position Covariance
% The robot and obstacle are simulated with process noise.
% An EKF is applied to estimate their states.
% After each update, the position covariance (for [px;py]) is extracted.

clear; clc; close all;

%% Simulation Parameters
dt = 0.05;            % time step [s]
T = 10;              % total simulation time [s]
N = T/dt;            % number of simulation steps

%% ---------------- Robot Setup ----------------
% Robot state: x_r = [px; py; psi]
x_r_true = [0; 0; 0];       % true initial state (may be inaccurate)
x_r_est  = [0; 0; 0];       % initial EKF estimate
P_r      = diag([0.1, 0.1, 0.1]);  % initial error covariance for robot

% Process noise covariance for robot dynamics
Q_r = diag([0.1, 0.1, 0.05]);

% Measurement noise covariance for robot (position measurement)
R_r = diag([0.6, 0.6].^2);

% Define a constant control input for the robot over time.
% Control: [vx; vy; omega]
u_r = repmat([1; 0; 0.05], 1, N);

% Storage for robot estimates and covariances
x_r_est_history = zeros(3, N);
P_r_history = zeros(3,3,N);

%% ---------------- Obstacle Setup ----------------
% Obstacle state: x_o = [px; py; vx; vy]
x_o_true = [5; 5; 0.5; 0];      % true initial state
x_o_est  = [5; 5; 0.5; 0];      % initial EKF estimate
P_o      = diag([1, 1, 0.1, 0.1]);  % initial error covariance for obstacle

% Process noise covariance for obstacle dynamics
Q_o = diag([0.1, 0.1, 0.05, 0.05]);

% Measurement noise covariance for obstacle (position measurement)
R_o = diag([0.3, 0.3]);

% Storage for obstacle estimates and covariances
x_o_est_history = zeros(4, N);
P_o_history = zeros(4,4,N);

%% ---------------- Simulation Loop ----------------
for k = 1:N
    %% --- True State Propagation with Process Noise ---
    % ---- Robot Propagation ----
    u = u_r(:, k);
    % Propagate true robot state using the discrete-time model and add process noise.
    x_r_true = [ x_r_true(1) + u(1)*dt;
                 x_r_true(2) + u(2)*dt;
                 x_r_true(3) + u(3)*dt ] + mvnrnd(zeros(3,1), Q_r)';
             
    % ---- Obstacle Propagation ----
    % Propagate true obstacle state using constant velocity model.
    x_o_true = [ x_o_true(1) + x_o_true(3)*dt;
                 x_o_true(2) + x_o_true(4)*dt;
                 x_o_true(3);
                 x_o_true(4) ] + mvnrnd(zeros(4,1), Q_o)';
    
    %% --- Measurements (Position Only) ---
    % For the robot: measure [px; py] with measurement noise.
    z_r = x_r_true(1:2) + mvnrnd(zeros(2,1), R_r)';
    % For the obstacle: measure [px; py] with measurement noise.
    z_o = x_o_true(1:2) + mvnrnd(zeros(2,1), R_o)';
    
    %% --- EKF for Robot ---
    % Prediction Step:
    % Predicted state using dynamics f(x,u) = [px+vx*dt; py+vy*dt; psi+omega*dt]
    x_r_pred = [ x_r_est(1) + u(1)*dt;
                 x_r_est(2) + u(2)*dt;
                 x_r_est(3) + u(3)*dt ];
    % Jacobian F for robot (here, identity because dynamics are linear in this simple model)
    F_r = eye(3);
    % Predicted covariance:
    P_r_pred = F_r * P_r * F_r' + Q_r;
    
    % Update Step:
    % Measurement model: z = H_r * x_r, where H_r extracts the position.
    H_r = [1 0 0;
           0 1 0];
    % Innovation:
    y_r = z_r - H_r * x_r_pred;
    % Innovation covariance:
    S_r = H_r * P_r_pred * H_r' + R_r;
    % Kalman Gain:
    K_r = P_r_pred * H_r' / S_r;
    % Updated state estimate:
    x_r_est = x_r_pred + K_r * y_r;
    % Updated covariance:
    P_r = (eye(3) - K_r * H_r) * P_r_pred;
    
    % Save robot estimates and covariance.
    x_r_est_history(:, k) = x_r_est;
    P_r_history(:,:,k) = P_r;
    
    %% --- EKF for Obstacle ---
    % Prediction Step:
    % Obstacle dynamics: f(x) = [px+vx*dt; py+vy*dt; vx; vy]
    x_o_pred = [ x_o_est(1) + x_o_est(3)*dt;
                 x_o_est(2) + x_o_est(4)*dt;
                 x_o_est(3);
                 x_o_est(4) ];
    % Jacobian F for obstacle:
    F_o = [1, 0, dt, 0;
           0, 1, 0, dt;
           0, 0, 1, 0;
           0, 0, 0, 1];
    % Predicted covariance:
    P_o_pred = F_o * P_o * F_o' + Q_o;
    
    % Update Step:
    % Measurement model: z = H_o * x_o, with H_o extracting the position.
    H_o = [1 0 0 0;
           0 1 0 0];
    % Innovation:
    y_o = z_o - H_o * x_o_pred;
    % Innovation covariance:
    S_o = H_o * P_o_pred * H_o' + R_o;
    % Kalman Gain:
    K_o = P_o_pred * H_o' / S_o;
    % Updated state estimate:
    x_o_est = x_o_pred + K_o * y_o;
    % Updated covariance:
    P_o = (eye(4) - K_o * H_o) * P_o_pred;
    
    % Save obstacle estimates and covariance.
    x_o_est_history(:, k) = x_o_est;
    P_o_history(:,:,k) = P_o;
end

%% ---------------- Extract and Display Position Covariance ----------------
% For the robot, extract the position covariance from P_r (submatrix corresponding to [px;py])
P_r_pos = P_r_history(1:2, 1:2, end);
fprintf('Robot Position Covariance Matrix at final time step:\n');
disp(P_r_pos);

% For the obstacle, extract the position covariance from P_o (submatrix corresponding to [px;py])
P_o_pos = P_o_history(1:2, 1:2, end);
fprintf('Obstacle Position Covariance Matrix at final time step:\n');
disp(P_o_pos);

%% ---------------- Plot Estimated Positions ----------------
figure;
plot(x_r_est_history(1,:), x_r_est_history(2,:), 'b.-','LineWidth',1.5); hold on;
plot(x_o_est_history(1,:), x_o_est_history(2,:), 'r.-','LineWidth',1.5);
legend('Robot Estimate','Obstacle Estimate');
xlabel('X Position [m]'); ylabel('Y Position [m]');
title('Estimated Positions from EKF');
grid on;
