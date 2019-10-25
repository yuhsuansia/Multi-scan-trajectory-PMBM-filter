clear;clc
dbstop if error

% Set simulation parameters
numtruth = 4; % number of targets
% pmi = 1; % simulation case 1 or 2
pmi = 2;
% covariance used for mid-point initialisation
Pmid = 1e-6*eye(4);
% detection probability
Pd = 0.9;
% clutter rate
lfai = 10;

% Algorithm setup
slideWindow = 3;
L = 1;

% Generate Ground Truth
[model,measlog,xlog,X] = gentruth(slideWindow,L,Pd,lfai,numtruth,Pmid,pmi);

% time steps in total
K = length(xlog);

% GOSPA metric
gospa_vals = zeros(K,4);
gospa_c = 10;

% Trajectory metric
traMetric = zeros(K,5);

% Simulation time
runningTimeperMCandTimeStep = zeros(K,1);

% Extract ground truth and re-construct for use of trajectory metric
Xtrue = cell(K,1);
for i = 1:K
    idx = i >= X.tVec;
    Xtrue{i}.tVec = X.tVec(idx);
    Xtrue{i}.iVec = min(X.iVec(idx),i-X.tVec(idx)+1);
    Xtrue{i}.xState = X.xState(:,1:max(Xtrue{i}.tVec+Xtrue{i}.iVec)-1,idx);
end

% Single trajectory hypothesis structure
%     r = existence probability;
%     xpre = predicted trajectory (target states across time);
%     x = updated trajectory (target states across time);
%     Ppre = predicted covariance;
%     P = updated covariance
%     l : trajectory (measurement indices across time)
%     c : the cost of each single target hypothesis
%     a : track index
%     beta: begin time of trajectory, absorbtion is used
%     epsilon: end time of trajectory, a vector
%     w = trajectory weight, a vector, w(i) specifies the weight of
%     trajectory with end time epsilon(i)

% Initialisation
trajectoryMBM = [];
unknownPPP = [];
trajectorySmoothEst = [];
xest = cell(K,1);
trajectory = cell(K,1);

% Loop through time
fprintf('Time step: ');
for t = 1:K
    
    fprintf('%g ', t);
    
    tstart = tic;
    
    % Predict all single target hypotheses of previous scan
    [trajectoryMBM,unknownPPP] = predictStep(trajectoryMBM,unknownPPP,model,t);
    
    % Update all predicted single target hypotheses of previous scan
    [unknownPPP,trajectoryUpdMBM,trajectoryNewMBM] = ...
        updateStep(unknownPPP,trajectoryMBM,model,measlog{t},t);
    
    % Multi-scan data association
    [trajectoryEst,trajectoryMBM] = dataAssoc(trajectoryUpdMBM,trajectoryNewMBM,model);
    
    if model.L ~= 0
        % Backforward Filtering for most likely global hypothesis
        trajectorySmoothEst = backforwardFiltering(trajectoryEst,trajectorySmoothEst,model);
        
        % Use part of smoothed trajectories to replace unsmoothed ones
        trajectoryMBM = approxSmoothedTrajectory(trajectoryMBM,trajectorySmoothEst,model);
    else
        trajectorySmoothEst = trajectoryEst;
    end
    
    % Target state extraction
    [xest{t},trajectory{t},Xest] = stateExtract(trajectorySmoothEst);
    
    % Record cycling time
    runningTimeperMCandTimeStep(t,1) = toc(tstart);
    
    % (Filtering) performance evaluation using GOSPA metric
    [d_gospa, ~, decomposed_cost] = GOSPA(get_comps(xlog{t},[1 3]), get_comps(xest{t},[1 3]), 2, gospa_c, 2);
    gospa_vals(t,:) = sqrt([d_gospa^2 decomposed_cost.localisation decomposed_cost.missed decomposed_cost.false]);
    
    % (Filtering or smoothing) performance evaluation using trajectory
    % metric, normalised by time step
    [dxy, ~, loc_cost, miss_cost, fa_cost, switch_cost] = LPTrajMetric_sparse(Xtrue{t}, Xest, gospa_c, 2, 2);
    traMetric(t,:) = sqrt([dxy^2/t sum(loc_cost)/t sum(miss_cost)/t sum(fa_cost)/t sum(switch_cost)/t]);
    
end

fprintf('\nAverage GOSPA: %g\nAverage Trajectory Metric: %g\n', mean(gospa_vals(:,1)), mean(traMetric(:,1)));
