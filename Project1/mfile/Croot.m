function F = Croot(matrix,n,M,C)
% function of solving C root when unknown C in a ML-CFAR.
% Copyright Xiangyu Ding, ECE in SSOE, University of Pittsburgh 2016
% Reference R. Ravid and N. Levanon, ¡°Maximum-likelihood CFAR for Weibull
% background,¡± Radar and Signal Processing, IEE Proceedings F, vol. 139,
% pp. 256-264, Jun. 1992.
v = matrix(n,:); % choose a specific pulse in raw data matrix
v3 = repmat(v,1,3); % cycle the pulse three times
col = numel(v); % compute the column number
for x = 1:col % compute every cell in the pulse
    s1(x) = 0; % 1s sum in the eqn.29
    s2(x) = 0; % 2nd sum in the eqn.29
    s3(x) = 0; % 3rd sum in the eqn.29
    for i = x+col-M/2:x+col+M/2 % compute the middle array
        s1(x) = s1(x) + (v3(1,i).^C(x))*log(v3(1,i)); % compute the 1st sum
        s2(x) = s2(x) + (v3(1,i).^C(x)); % compute the 2nd sum
        s3(x) = s3(x) + log(v3(1,i)); % compute the 3rd sum
    end
    s1(x) = s1(x) - (v3(1,x).^C(x))*log(v3(1,x)); % minus the CUT cell
    s2(x) = s2(x) - (v3(1,x).^C(x));
    s3(x) = s3(x) - log(v3(1,x));
    F(x) = (s1(x)/s2(x)-s3(x)/M)-1/C(x); % build the function for parameter
    % C
end