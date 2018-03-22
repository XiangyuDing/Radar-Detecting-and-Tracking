function [x_Detect,T] = Cfarx(matrix,Pfa,M)
% function of ML-CFAR detector in known C=2 and total data.
% Copyright Xiangyu Ding, ECE in SSOE, University of Pittsburgh 2016
% Reference R. Ravid and N. Levanon, ¡°Maximum-likelihood CFAR for Weibull
% background,¡± Radar and Signal Processing, IEE Proceedings F, vol. 139,
% pp. 256-264, Jun. 1992.
row = size(matrix,1); % compute the number of pulse
col = size(matrix,2); % compute the number of column
for n = 1:row
    v = matrix(n,:); % choose a specific pulse in raw data matrix
    v3 = repmat(v,1,3); % cycle the pulse three times
    for x = 1:col % compute every cell in the pulse
        y = x + col; % change the cell into the middle array
        sum = 0; % define sum
        for i = y-M/2:y+M/2 % calculate the sum
            sum = sum + v3(1,i)^2;
        end
        sum = sum - v3(1,y)^2; % minus the CUT cell
        B = (sum/M)^(1/2); % calculate B
        alpha = (((Pfa^(-1/M))-1)*M)^(1/2); % calculate alpha
        t = abs(alpha * B); % calculate threshold in cell x
        if v3(1,x)>t  % input exceeds threshold
            x_detect = 1; % target detected
        else
            x_detect = 0; % no target detected
        end
        x_Detect(n,x) = x_detect; % save the total detection results of data
        T(n,x) = t; % save the total threshold of data
    end
end