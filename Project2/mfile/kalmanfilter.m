function X = kalmanfilter(T,k,x0,Zm,option)
% function of Kalman filter to track the target motion
% Input -- sample period T; sample times k; initial condition x0;
% measured location Zm;
% reference model type option(1 is constant velocity model and 2 is
% coordinated turn model).
% Output -- target motion condition X.
% Copyright Xiangyu Ding, ECE in SSOE, University of Pittsburgh Dec,2016
% Reference --
% https://www.youtube.com/watch?v=2-lu3GNbXM8&t=861s

x = x0; % initial condition
n = pi / k; % turn rate omega

X = zeros(4,4*k); % initial total condition matrix
P = 0.001*eye(4); % initial state covariance matrix
Fcv = [1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1]; % constant velocity model
Fct = [1,sin(n*T)/n,0,-(1-cos(n*T))/n;0,cos(n*T),0,-sin(n*T);0,(1-cos(n*T))/n,1,sin(n*T)/n;0,sin(n*T),0,cos(n*T)]; % coordinated turn model
Q = 0.001*eye(4); % motion model covariance matrix
H = [1,0,0,0;0,0,1,0]; % measurement model
R = [0.1,0;0,0.1]; % measurement noise covariance matrix

switch option % select the reference model
    case 1 % reference model is constant velocity model
        F = Fcv;
        
    case 2 % reference model is coordinated turn model
        F = Fct;
end

for i = 1:4*k % filter the total course
    x_ = F*x; % estimate location
    P_ = F*P*F'+Q; % estimate the state covariance matrix
    S = H*P_*H'+R; % innovation covariance matrix
    K = P_*H'/S; % Kalman operator equation
    x = x_+K*(Zm(:,i)-H*x_); % update the estimating location
    P = (eye(4)-K*H)*P_; % update the state covariance matrix
    X(:,i) = x; % save the current location into the total condition matrix
end