%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Written by : Casper Versteeg & Alex Lin
%   Date: 3/30/2018
%   University of Georgia SSRL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%INITIALIZING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = [0; 0; 0; 1];                   % Initial attitude quaternion
w = [0.1; 0.1; 0.1];                % Initial satellite spin
q_set = [0; 0; 0; 1];               % Attitude setpoint

epoch = [2015; 1; 1; 0; 0; 0];      % Epoch time for Julian Date
N = 4;                              % Number of orbits to model

t = dlmread('time.42');             % Vector containing time since epoch
PosN = dlmread('PosN.42');          % Position vector in inertial coordinates at t
VelN = dlmread('VelN.42');          % Velocity vector in inertial coordinates at t
dt = t(2)-t(1);                     % Timestep in 't' vector
refreshRate = 1;                    % ADCS frequency
% Filter time, position and velocity to ADCS refresh rate
t = 1:refreshRate:N*t(end);
PosN = PosN(1:refreshRate/dt:end,:);
VelN = VelN(1:refreshRate/dt:end,:);

J = [3.98e-2 1.96e-4 2.1e-5;
    1.96e-4 7.54e-3 -9.98e-5; 
    2.1e-5 -9.98e-5 3.96e-2];       % Satellite inertia matrix
Jw = diag(ones(1,3)*3.534e-3/(7500*2*pi/60)); % Reaction wheel inertia
Ai = eye(3);                        % Reaction wheel layout matrix
wrw = [0; 0; 0];                    % Initial reactio wheel spin
TAM = eye(3);                       % Magnetometer layout matrix
FSS1 = [-1 0 0; 0 1 0; 0 0 -1]';    % Fine sun sensor layout matrix (North)
FSS2 = [0 0 -1 ;0 1 0 ;-1 0 0]';    % Fine sun sensor layout matrix (Zenith)
QFSS = (1)^2;                        % Fine sun sensor covariance
QTAM = (1)^2;                        % Magnetometer covariance
Qgyr = (1)^2;                        % Gyroscope covariance
Q = diag([QFSS, QFSS, QTAM, Qgyr]); % Covariance matrix

B_meas = zeros(3, length(t));       % Initialize measurements
sF1 = zeros(3, length(t)); 
sF2 = zeros(3, length(t)); 
w_meas = zeros(3, length(t)); 
q_B = zeros(4, length(t));          % Initialize attitude estimates
q_F1 = zeros(4, length(t)); 
q_F2 = zeros(4, length(t)); 
q_w_dot = zeros(4, length(t)); 
q_w = zeros(4, length(t));
%%%%%%%%%%%%%%%%%%%%%%%%END INITIALIZING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start Kalman Filter
EKF = extendedKalmanFilter(@stateTransFcn, @measurementStateFcn, q);
EKF.MeasurementNoise = Q;

% Build attitude and space propagator model
for i = 1:length(t)
    % Create renormalized index that accounts for orbit propagator past
    % the end of 'PosN'
    I = i - length(PosN)*floor(t(i)/length(PosN));
    % Compute environment variables: magnetic field, sun vector, 
    % aerodynamic disturbance torque, radiation pressure torque and 
    % gravity gradient torque.
    [~,A] = quaternionAttitudeMatrix(q);
    %B = [];
    [Tsol, ~, s] = solarRadPress(PosN(I,:), A, t(i), 'cart', epoch);
    Taero = aerodynDrag(VelN(I,:), A);
    %Tgrav = [];
    
    % Find sensor measurements from Fine sun sensors, magnetometer and
    % gyroscope
    sF1(:,i) = fineSunSensor(s, A, FSS1);
    sF2(:,i) = fineSunSensor(s, A, FSS2);
    B_meas(:,i) = threeAxisMagnetometer(B, A, TAM);
    w_meas(:,i) = gyroscope(w);
    % Obtain attitude estimations from these measurements
    q_F1(:,i) = quaternionAttitudeMatrix(attitudeMatrix(sF1(:,i), s));
    q_F2(:,i) = quaternionAttitudeMatrix(attitudeMatrix(sF2(:,i), s));
    q_B(:,i) = quaternionAttitudeMatris(attitudeMatrix(B_meas(:,i), B));
    q_w_dot(:,i) = quaternionRate(q, w_meas);
    q_w(:,i) = q_w_dot(:,i)*refreshRate+q_B(:,i);
    
end


% Define attitude state transition function:
function [q, w] = stateTransFcn(q, w)
    J = [3.98e-2 1.96e-4 2.1e-5;
        1.96e-4 7.54e-3 -9.98e-5; 
        2.1e-5 -9.98e-5 3.96e-2];
    Jw = diag(ones(1,3)*3.534e-3/(7500*2*pi/60));
    Ai = eye(3); dt = 1;
    J_hat = J-(Ai*Jw)*Ai';
    w_dot = J_hat\(-crossMatrix(w)*(J*w+Ai*Jw*h)+Ai*Tc+Tm+Td);
    w = w+w_dot*dt;
    q_dot = quaternionRate(q, w);
    q = q+q_dot*dt;
end
% Define measurement state function
function [q_F1, q_F2, q_B, q_w] = measurementStateFcn(q_F1, q_F2, q_B, q_w)
    
end
