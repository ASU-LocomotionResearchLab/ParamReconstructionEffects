function [out] = investigateFNN(x, dataLeng, taus, alt)
%investigateFNN - Investigating effect of different values of tau on False
%Nearest Neighbors (FNN) output: embedding dimension (mDim). This function
%will determine multiple mDim for each of the tau inputs, looking at
%different threshold values for consistency.
% ========================================================================
% INPUTS:
%   x    - Nx1 data
%   taus - Nx1 array holding different threshold values for tau
%           Ex. order for LES Study:
%          (1) 1st min tau (2) tau < 1/e (3) tau < 0.1 (4) tau < 0.05
%   alt  = used for alternate testing schema, enter 1 for this to be run
%
% OUTPUTS:
%   out - struct that contains the following for each tau in the taus vector:
%         fnn, tau, mDim_min, mDim_minp, mDim_10p, mDim_5p, mDim_1p
%
% NOTES:
% 1. Threshold values for determining mDim are:
%    (1) minimum (2) mDim < 10% (3) mDim < 5% (4) mDim < 1%
%    -- These values are ideal for studying gait and therefore might need
%       to be changed if studying a different system
% 2. If any value of mDim is set to 0, then no values crossed the
%    threshold to be counted
%
%
% Created By: Victoria Smith Hussain
% Created On: 1/29/2019
% ========================================================================
% If alt is not used, alt = 0
if ~exist('alt', 'var')
    alt = 0;
end

% Check data and make a raw vector of scalar data
if min(size(x)) ~= 1
    error('Data input should be a vector of scalar data points.')
else
    x=x(:)'; % just the way we need it
end
% Initializing variables
maxDim = 8;

if alt == 0
    % Check taus and make a raw vector of scalar data
    if min(size(taus)) ~= 1
        error('Taus should be a Nx1 vector of scalar data points.')
    else
        taus=taus(:)'; % just the way we need it
    end
    
    % Loop
    tausNames = {'tau_min'; 'tau_e'; 'tau_10p'; 'tau_5p'};
    fprintf(' (')
    for ii = 1:length(taus)
        % Determine the max data length possible for FNN
        maxDataLeng = dataLeng - maxDim*taus(ii);
        
        fprintf(' %d ',ii)
        if taus(ii) == 0
            mDim_10p = 0; mDim_5p = 0;
            mDim_1p = 0;  fnn = []; deltaFNN = []; mDim_d = 0;
        else
            % Calculate FNN
            % Defaults in fnns: R_tol = 15; A_tol = 4;
            [fnn] = fnns(x,taus(ii),maxDim,maxDataLeng);
            deltaFNN = -fnn(2:end) + fnn(1:end-1);
            [~, ~, mDim_10p, mDim_5p, mDim_1p, mDim_d] = findEmbeddingDims(fnn, deltaFNN);
        end
        
        
        
        % Save data before next iteration & for the final output
        out(ii).fnn        = fnn;
        out(ii).deltaFNN   = deltaFNN;
        out(ii).tau        = taus(ii);
%         out(ii).mDim_min   = mDim_min;
%         out(ii).mDim_minp  = mDim_minp;
        out(ii).tType      = string(tausNames(ii));
        out(ii).mDim_d     = mDim_d;
        out(ii).mDim_10p   = mDim_10p;
        out(ii).mDim_5p    = mDim_5p;
        out(ii).mDim_1p    = mDim_1p;
       
    end
    fprintf(') ')
    
elseif alt == 1
    if length(taus) == 2
        tausNames = {'ideal'; 'arb'};
    elseif length(taus) == 3
        tausNames = {'sets';'ideal'; 'arb'};
    else
        disp('No titles for taus: labeled 1-2-3-4')
        tausNames = {'1'; '2'; '3'; '4'};
    end
    
    % Duplicates?
    [~,unqInd] = unique(taus);

    
    fprintf(' (')
    for ii = 1:length(taus)
        if sum(ii == unqInd)
            
            fprintf(' %d ',ii)
            % Determine the max data length possible for FNN
            maxDataLeng = dataLeng - maxDim*taus(ii);
            
            if taus(ii) == 0
                mDim_min = 0; mDim_minp = nan; mDim_10p = 0; mDim_5p = 0;
                mDim_1p = 0;  fnn = []; mDim_d = 0;
            else
                % Calculate FNN
                % Defaults in fnns: R_tol = 15; A_tol = 2;
                [fnn] = fnns(x,taus(ii),maxDim,maxDataLeng);
                deltaFNN = -fnn(2:end) + fnn(1:end-1);
                [~, ~, mDim_10p, mDim_5p, mDim_1p, mDim_d] = findEmbeddingDims(fnn,deltaFNN);
            end
            
            % Save data before next iteration & for the final output
            out(ii).fnn        = fnn;
            out(ii).deltaFNN   = deltaFNN;
            out(ii).tau        = taus(ii);
            out(ii).tType      = string(tausNames(ii));
            out(ii).mDim_delta = mDim_d;
            %         out(ii).mDim_min   = mDim_min;
            %         out(ii).mDim_minp  = mDim_minp;
            out(ii).mDim_10p   = mDim_10p;
            out(ii).mDim_5p    = mDim_5p;
            out(ii).mDim_1p    = mDim_1p;
        
        else
            fprintf(' copy ')
            out(ii).fnn        = out(ii-1).fnn;
            out(ii).deltaFNN   = out(ii-1).deltaFNN;
            out(ii).tau        = taus(ii);
            out(ii).tType      = string(tausNames(ii));
            out(ii).mDim_delta = out(ii-1).mDim_delta;
            %         out(ii).mDim_min   = mDim_min;
            %         out(ii).mDim_minp  = mDim_minp;
            out(ii).mDim_10p   = out(ii-1).mDim_10p;
            out(ii).mDim_5p    = out(ii-1).mDim_5p;
            out(ii).mDim_1p    = out(ii-1).mDim_1p;
            
        end
        
    end
    fprintf(') ')
    
    
end
end

function [mDim_min, mDim_per, mDim_10p, mDim_5p, mDim_1p, mDim_d] = findEmbeddingDims(fnn, dfnn)


%------------ Finding the Minimum
[pks,ipks] = findpeaks(-fnn);
if ~isempty(ipks)
    mDim_min = ipks(1);
    mDim_per = -pks(1);
else
    [mDim_per, mDim_min] = min(fnn);
end

%------------ Finding first point < 10%
temp = fnn < 0.1;
if sum(temp) > 0
    firstPt = find(temp); firstPt = firstPt(1);
    mDim_10p = firstPt;
    clearvars firstPt temp
else
    mDim_10p = 0;
    clearvars temp
end

%------------ Finding first point < 5%
temp = fnn < 0.05;
if sum(temp) > 0
    firstPt = find(temp); firstPt = firstPt(1);
    mDim_5p = firstPt;
    clearvars firstPt temp
else
    mDim_5p = 0;
    clearvars temp
end

%------------ Finding first point < 1%
temp = fnn < 0.01;
if sum(temp) > 0
    firstPt = find(temp); firstPt = firstPt(1);
    mDim_1p = firstPt;
    clearvars firstPt temp
else
    mDim_1p = 0;
    clearvars temp
end

%------------ Finding first point < 5% in delta FNN k
% We want to find when the change between dimensions is less than 5%
temp = dfnn < 0.05 & dfnn > 0;
if sum(temp) > 0
    firstPt = find(temp); myPoint = firstPt(1);
    if myPoint == 1
        myPoint = firstPt(2);
    end
    mDim_d = myPoint; 
    clearvars firstPt temp
else
    mDim_d = 0;
    clearvars temp
end

end
