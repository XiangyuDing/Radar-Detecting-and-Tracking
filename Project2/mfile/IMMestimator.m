function [X,U] = IMMestimator(T,k,x0,Zm)
% function of The Interacting Multiple Model Estimator Processing
% Input -- sample period T; sample times k; initial condition x0;
% measured location Zm;
% Output -- target motion condition X; mode probability U.
% Copyright Xiangyu Ding, ECE in SSOE, University of Pittsburgh Dec,2016
% Reference --
% Radar Tracking with an Interacting Multiple Model and
% Probabilistic Data Association Filter for Civil
% Aviation Applications,Shau-Shiun Jan * and Yu-Chun Kao

% initial parameter
x1 = x0; % initial condition
x2 = x0;
P1 = 0.001*eye(4); % initial state covariance matrix
P2 = 0.001*eye(4);
n = pi / k; % turn rate omega

Fcv = [1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1]; % constant velocity model
Fct = [1,sin(n*T)/n,0,-(1-cos(n*T))/n;0,cos(n*T),0,-sin(n*T);0,(1-cos(n*T))/n,1,sin(n*T)/n;0,sin(n*T),0,cos(n*T)]; % coordinated turn model
Q = 0.001*eye(4); % motion model covariance matrix
H = [1,0,0,0;0,0,1,0]; % measurement model
R = [0.1,0;0,0.1]; % measurement noise covariance matrix

u1 = 0.5; % initial probability of mode1
u2 = 0.5; % initial probability of mode2
p = [0.95,0.05;0.05,0.95]; % Markov chain matrix
X = zeros(4,4*k); % initial condition output matrix
U = zeros(2,4*k); % initial mode probability matrix

for i =1:4*k
    % calculation of the mixing probabilities
    c1 = p(1,1)*u1+p(2,1)*u2;
    c2 = p(1,2)*u1+p(2,2)*u2;
    u11 = p(1,1)*u1/c1;
    u12 = p(1,2)*u1/c2;
    u21 = p(2,1)*u2/c1;
    u22 = p(2,2)*u2/c2;
    % mixing
    x01 = x1*u11+x2*u21;
    x02 = x1*u12+x2*u22;
    P01 = u11*(P1+(x1-x01)*(x1-x01)')+u21*(P2+(x2-x01)*(x2-x01)');
    P02 = u12*(P1+(x1-x02)*(x1-x02)')+u22*(P2+(x2-x02)*(x2-x02)');
    % Kalman filter implementation
    x1_ = Fcv*x01;
    P1_ = Fcv*P01*Fcv'+Q;
    K1 = P1_*H'/(H*P1_*H'+R);
    x1 = x1_+K1*(Zm(:,i)-H*x1_);
    P1 = (eye(4)-K1*H)*P1_;
    x2_ = Fct*x02;
    P2_ = Fct*P02*Fct'+Q;
    K2 = P2_*H'/(H*P2_*H'+R);
    x2 = x2_+K2*(Zm(:,i)-H*x2_);
    P2 = (eye(4)-K2*H)*P2_;
    % mode-matched filtering
    S1 = H*P1_*H'+R;
    y1 = Zm(:,i)-H*x1_;
    d1 = y1'*S1^-1*y1;
    A1 = exp(-0.5*d1)/sqrt(det(2*pi*S1));
    S2 = H*P2_*H'+R;
    y2 = Zm(:,i)-H*x2_;
    d2 = y2'*S2^-1*y2;
    A2 = exp(-0.5*d2)/sqrt(det(2*pi*S2));
    % mode probability update
    c = A1*c1+A2*c2;
    u1 = A1/c*c1;
    u2 = A2/c*c2;
    % estimate and covariance combination
    x = x1*u1+x2*u2;
    P = u1*(P1+(x1-x)*(x1-x)')+u2*(P2+(x2-x)*(x2-x)');
    % save the estimation and mode probability
    X(:,i) = x;
    U(:,i) = [u1;u2];
end