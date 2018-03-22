function [x_Detect,T] = cfar(matrix,n,Pfa,M)
% function of ML-CFAR detector in known C=2 and specific pulse. 
% Copyright Xiangyu Ding, ECE in SSOE, University of Pittsburgh 2016
% Reference R. Ravid and N. Levanon, ¡°Maximum-likelihood CFAR for Weibull 
% background,¡± Radar and Signal Processing, IEE Proceedings F, vol. 139,
% pp. 256-264, Jun. 1992.
v = matrix(n,:); % choose a specific pulse in raw data matrix
v3 = repmat(v,1,3); % cycle the pulse three times
col = numel(v); % compute the column number
for x = 1:col % compute every cell in the pulse
    x = x + col; % change the cell into the middle array
    sum = 0; % define sum
    for i = x-M/2:x+M/2 % calculate the sum
        sum = sum + v3(1,i)^2;
    end
    sum = sum - v3(1,x)^2; % minus the CUT cell
    B = (sum/M)^(1/2); % calculate B
    alpha = (((Pfa^(-1/M))-1)*M)^(1/2); % calculate alpha
    t = abs(alpha * B); % calculate threshold in cell x
    if v3(1,x)>t  % input exceeds threshold
        x_detect = v3(1,x); % target detected
    else
        x_detect = nan; % no target detected
    end
    x_Detect(x) = x_detect; % save the total detection results of pulse
    T(x) = t; % save the total threshold of pulse
end