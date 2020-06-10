function [taus, pDist] = investigateAMI(x)
%investigateAMI - Investigating different thresholds for defining tau using
%Average Mutual Information to be used in False Nearest Neighbors (FNN) and
%subsequent NLD analysis (e.g. Lyapunov Exponents)
% =========================================================================
% INPUTS:
%   data - Nx1 array
%
% OUTPUTS:
%   taus  - array holding different threshold values for tau in the order:
%           (1) 1st min tau (2) tau < 1/e (3) tau < 0.1 (4) tau < 0.05
%   pDist -
% 
% NOTES:
%   1. If any value of tau is set to 0, then no values crossed the
%   threshold to be counted
%
%
% Created By: Victoria Smith Hussain - vasmith5@asu.edu
% Created On: 1/29/2019
% =========================================================================

% check data and make a raw vector of scalar data
if min(size(x)) ~= 1
    error('Input should be a vector of scalar data points.')
else
    x=x(:)'; % just the way we need it
end

% Calculate AMI given x
pDist = []; %clear all contents of v
for i=0:100
    [pDist(i+1),lag] = ami(x,x,i);
end

% First Minimum tau
[~,ipks] = findpeaks(-pDist);
if isempty(ipks)
    [~, tau_min] = min(pDist);
else
    tau_min = ipks(1);
end

% Find tau when < 1/e (~0.3679)
temp = pDist < 1/exp(1);
if sum(temp) > 0
    firstPt = find(temp);
    tau_e = firstPt(1);
    clearvars temp firstPt
else 
    tau_e = 0;
    clearvars temp
end

% Find tau when < .1
temp = pDist < .1;
if sum(temp) > 0
    firstPt = find(temp);
    tau_10p = firstPt(1);
    clearvars temp firstPt
else 
    tau_10p = 0;
    clearvars temp
end

% Find tau when < .05
temp = pDist < 0.05;
if sum(temp) > 0
    firstPt = find(temp);
    tau_5p = firstPt(1);
    clearvars temp firstPt
else 
    tau_5p = 0;
    clearvars temp
end

% Final outputs
taus = [tau_min tau_e tau_10p tau_5p];
pDist = pDist';

end

