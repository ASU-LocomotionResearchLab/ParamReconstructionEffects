%=========================================================================%
%                   LES_Tau_Dim_DataLeng_Inv.m
%
% Need to investigate tables for storing results at a certain level to it's
% not such a nested Results struct
% ------------------------- APDM Orientation ---------------------------- %
%           |----------------|
%           | O     ^ x      |          x (1) - Vertical
%           |       |        |          y (2) - ML Direction
%           |   <-- o        |          z (3) - AP Direction
%           |     y          |
%           |----------------|
%       Z is straight through device
% ----------------------------------------------------------------------- %
%
%   Created By: Victoria Smith Hussain
%   Created On: 1/29/2019
%   Updated On: 3/04/2019
%
%   Based on LES_Investigation_Tau_mDim_DataLeng & LES_Inv_mDim_DataLeng
%=========================================================================%

clear; close all; clc;

% addpath(genpath('D:\Matlab_Workspace\Dissertation\LES'))

%Load subject data
load('LES_LTW.mat');
fs = 128;
% Determine the # of subjects to be processed
nSubj = numel(s);

% Initiation Vars for how to evaluate the role of data length in AMI and FNN
%   gcSet is used for raw & gcNorm # of Strides
%   dpSet is for # of points.
gcSet = [30 50 100 150 200 300 500 700 900 1100];

% Based on 15s,30s,1m,2m,3m,5m,10m,15m,20m of data when fs=128
dpSet = [2000 4000 7700 15000 23000 38000 77000 115000 153000];%

% Naming Conventions
saveFilesNames = {'R_LES_InvTauDim_gc';'R_LES_InvTauDim_gcNorm';'R_LES_InvTauDim_dp'};
dirNames = {'VT'; 'AP'; 'ML'};
typeName = {'gc'; 'gcNorm'; 'dp'};
arbTau = 10;

%% =======================================================================%
%                                TESTING AMI & FNN
% - Finding how tau changes wrt gait data length
%   - Based on the # of strides, # of normalized strides, and # of data
%      points (w/o respect to strides)
%   - Collecting information about the minimum tau, tau < 1/e, tau < 10%,
%     and tau < 5%
% - Finding how mDim changes wrt gait data length
%   - Based on the # of strides, # of normalized strides, and # of data
%      points (w/o respect to strides)
%   - Calculating mDim based on each signal's tau found using AMI, and then
%       using an arbitrary tau = 10
% ========================================================================%
type   = 2; %1-gc; 2-gcNorm; 3-dp
stSubj = 1;
enSubj = 1;

idealGC = 6;
idealDP = 4; 

finalSaveFileName = strcat(string(saveFilesNames(type)), '_s', string(stSubj),'-',string(enSubj), '.mat');

switch string(typeName(type))
    % ========================================================================
    %
    %
    % ========================================================================
    case {'gc', 'gcNorm'}
        
        for subj = stSubj:enSubj
            fprintf('Subj:%d',subj);
            
            for dir = 1:3 % VT;AP;ML (in order)
                fprintf('\n dir:%d\t AMI gc:',dir);
                HC   = s(subj).t(type).HC.data;
                
                
                %==========================================================
                %                     AMI INVESTIGATION
                %==========================================================
                for gc = 1:7 %length(gcSet)
                    tic
                    fprintf(' %d',gc);
                    % -------- Saving Data --------
                    % Load file if already exists to append new data, if not create the
                    % savefile with first data set
                    if(exist(finalSaveFileName, 'file'))
                        load(finalSaveFileName);
                    end
                    
                    % Check if max numStrides or numPts is too large before
                    % running analysis
                    if gc > length(HC)
                        continue %go to the next iteration of gc
                    else
                        
                        numStride = gcSet(gc);
                        % Extract desired time series
                        if type == 1
                            data = s(subj).t(type).v(dir).data(1:HC(numStride*2)+1);
                        else
                            data = s(subj).t(type).v(dir).data(1:HC(numStride)+1);
                        end
                        % Calculate various possible taus
                        [taus, pDist] = investigateAMI(data);
                        
                        tauTime = toc;
                        % Organzie output struct
                        R(subj).name        = s(subj).name;
                        R(subj).type        = string(typeName(type));
                        R(subj).v(dir).name = string(dirNames(dir));
                        R(subj).v(dir).set(gc).name     = strcat(num2str(gcSet(gc)), string(typeName(type)));
                        R(subj).v(dir).set(gc).pDist    = pDist;
                        R(subj).v(dir).set(gc).tauMin   = taus(1);
                        R(subj).v(dir).set(gc).tTime    = tauTime;
                        
                        % -------- Saving Data --------
                        if subj > stSubj
                            save(finalSaveFileName, 'R', '-append')
                        else
                            save(finalSaveFileName, 'R')
                        end
                        % --------
                        
                        clearvars taus pDist data numStride
                    end
                end % END of GC LOOP - 1
                
                
                %==========================================================
                %                     FNN INVESTIGATION
                %==========================================================
                
                for gc = 1:7 %length(gcSet)
                    tic
                    fprintf('\n\t\t FNN gc: %d',gc');
                    % Check if max numStrides or numPts is too large before
                    % running analysis
                    if gc > length(HC)
                        continue %go to the next iteration of gc
                    else
                        % Extract desired time series
                        numStride = gcSet(gc);
                        data = s(subj).t(type).v(dir).data(1:HC(numStride*2));
                        % Define taus to be tested: 1) that set's tau,
                        % 2) ideal tau: tau @gc = 300; 3) arbTau defined at top
                        taus = [R(subj).v(dir).set(gc).tauMin, R(subj).v(dir).set(idealGC).tauMin, arbTau];
                        % Calculate various possible mDims based on taus
                        [outFNN] = investigateFNN(data, length(data), taus, 1);
                        
                        dimTime = toc;
                        % Organzie output struct
                        R(subj).v(dir).set(gc).dimData = outFNN;
                        R(subj).v(dir).set(gc).dTime   = dimTime;
                        
                        % -------- Saving Data --------
                        save(finalSaveFileName, 'R', '-append')
                        
                        clearvars taus outFNN data numStride
                    end
                end
                
            end % END OF DIR LOOP
            fprintf('\n');
        end % END OF SUBJ LOOP
        % ========================================================================
        %
        %
        % ========================================================================
    case 'dp'
        for subj = stSubj:enSubj
            fprintf('Subj:%d',subj);
            
            for dir = 1:3 % VT;AP;ML (in order)
                fprintf('\n dir:%d\t AMI dp:',dir);
                
                %==========================================================
                %                     AMI INVESTIGATION
                %==========================================================
                for dp = 1:7 %length(dpSet)
                    tic
                    fprintf(' %d',dp);
                    % -------- Saving Data --------
                    % Load file if already exists to append new data, if not create the
                    % savefile with first data set
                    if(exist(finalSaveFileName, 'file'))
                        load(finalSaveFileName);
                    end
                    
                    numPts = dpSet(dp);
                    maxLeng = length(s(subj).t(1).v(dir).data);
                    
                    % Check if numPts is too large before running analysis
                    if numPts > maxLeng
                        continue % Go to the next iteration of dp
                    else
                        % Extract desired time series
                        data = s(subj).t(1).v(dir).data(1:numPts);
                        
                        % Calculate various possible taus
                        [taus, pDist] = investigateAMI(data);
                        
                        tauTime = toc;
                        % Organzie output struct
                        R(subj).name        = s(subj).name;
                        R(subj).type        = string(typeName(type));
                        R(subj).v(dir).name = string(dirNames(dir));
                        R(subj).v(dir).set(dp).name   = strcat(num2str(dpSet(dp)), string(typeName(type)));
                        R(subj).v(dir).set(dp).pDist  = pDist;
                        R(subj).v(dir).set(dp).tauMin = taus(1);
                        R(subj).v(dir).set(dp).tTime  = tauTime;
                        clearvars data numPts taus pDist tauTime
                        
                        % -------- Saving Data --------
                        if subj > stSubj
                            save(finalSaveFileName, 'R', '-append')
                        else
                            save(finalSaveFileName, 'R')
                        end
                        % --------
                        
                    end
                    
                end
                

                
                %==========================================================
                %                     FNN INVESTIGATION
                %==========================================================
                           
                for dp = 1:7 %length(dpSet)
                    tic
                    fprintf('\n\t\t FNN dp: %d',dp');   
                    
                    numPts = dpSet(dp);
                    maxLeng = length(s(subj).t(1).v(dir).data);
                    
                    % Check if numPts is too large before running analysis
                    if numPts > maxLeng
                        continue % Go to the next iteration of dp
                    else
                        % Extract desired time series
                        data = s(subj).t(1).v(dir).data(1:numPts);
                        
                        % Define taus to be tested: 1) that set's tau,
                        % 2) ideal tau: tau @gc = 300; 3) arbTau defined at top
                        taus = [R(subj).v(dir).set(dp).tauMin, R(subj).v(dir).set(idealDP).tauMin, arbTau];
                        % Calculate various possible mDims based on taus
                        [outFNN] = investigateFNN(data, length(data), taus, 1);
                        dimTime = toc;
                        
                        % Organzie output struct
                        R(subj).v(dir).set(dp).dimData = outFNN;
                        R(subj).v(dir).set(dp).dTime   = dimTime;
                        
                        % -------- Saving Data --------
                        save(finalSaveFileName, 'R', '-append')
                        
                    end
                end
                
            end % END of DIR LOOP
            fprintf('\n');
        end % END OF SUBJ LOOP
        
end % END OF SWITCH-CASE
