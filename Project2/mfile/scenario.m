function [Zt,Zm] = scenario(T,k,x0)
% function of scenario to simulate the target motion
% input -- sample period T; sample times k; initial condition x0; 
% output -- target motion true location Zt and measured location Zm.
% Copyright Xiangyu Ding, ECE in SSOE, University of Pittsburgh Nov,2016
% Reference Lab Assignment 2

x = x0; % initial condition
n = pi / k; % turn rate omega

X = zeros(4,4*k); % initial total condition matrix
Fcv = [1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1]; % constant velocity model
Fct = [1,sin(n*T)/n,0,-(1-cos(n*T))/n;0,cos(n*T),0,-sin(n*T);0,(1-cos(n*T))/n,1,sin(n*T)/n;0,sin(n*T),0,cos(n*T)]; % coordinated turn model
N = [0.5*T^2,0;T,0;0,0.5*T^2;0,T]; % input noise matrix
H = [1,0,0,0;0,0,1,0]; % measurement model

% target motion process
for i = 1:1:k % first loop
    x = Fcv * x + N * normrnd(0,sqrt(0.001),[2,1]); % target motion with white Gaussian noise
    X(:,i) = x; % save the current location into the total condition matrix
end
for i = k+1:1:2*k  % second loop
    x = Fct * x + N * normrnd(0,sqrt(0.001),[2,1]);
    X(:,i) = x;
end
for i = 2*k+1:1:3*k % third loop
    x = Fcv * x + N * normrnd(0,sqrt(0.001),[2,1]);
    X(:,i) = x;
end
for i = 3*k+1:1:4*k % forth loop
    x = Fct * x + N * normrnd(0,sqrt(0.001),[2,1]);
    X(:,i) = x;
end

Zt = H * X;  % true course

Zm = Zt + normrnd(0,sqrt(0.1),[2,4*k]); % measured course(add white Gaussian noise)