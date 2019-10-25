function [trajectoryMBM,MB_birth] = predictStep(trajectoryMBM,model,t)
%PREDICT: PREDICT MULTI-BERNOULLI AND POISSON COMPONENTS

% Get multi-Bernoulli prediction parameters from model
F = model.F;
Q = model.Q;
Ps = model.Ps;

% Interpret length of inputs
nb = length(model.rb);      % number of bernoulli components in MB birth
n = length(trajectoryMBM);

% Implement prediction algorithm

% Predict existing tracks (single target hypotheses)
for i = 1:n
    % if the probability that the trajectory still exists at current time
    % step is very small, then there is no need to predict it
    if trajectoryMBM(i).r~=0 && trajectoryMBM(i).w(end) > model.threshold
        
        xpre = F*trajectoryMBM(i).x(:,end);
        Ppre = F*trajectoryMBM(i).P(:,:,end)*F'+Q;
        
        % trajectory case
        trajectoryMBM(i).x = [trajectoryMBM(i).x xpre];
        trajectoryMBM(i).P = cat(3,trajectoryMBM(i).P,Ppre);
        
        % store prediction
        trajectoryMBM(i).xpre = [trajectoryMBM(i).xpre xpre];
        trajectoryMBM(i).Ppre = cat(3,trajectoryMBM(i).Ppre,Ppre);
        
        % beta remains the same, only update epsilon
        trajectoryMBM(i).epsilon = [trajectoryMBM(i).epsilon trajectoryMBM(i).epsilon(end)+1];
        
        % update weight vector of death time
        trajectoryMBM(i).w = [trajectoryMBM(i).w(1:end-1) trajectoryMBM(i).w(end)*(1-Ps) trajectoryMBM(i).w(end)*Ps];
    end
end

% Incorporate birth intensity
MB_birth = repmat(struct,nb,1);
% Allocate memory
for i = 1:nb
    MB_birth(i).r = model.rb(i);
    MB_birth(i).x = model.xb(:,i);
    MB_birth(i).P = model.Pb(:,:,i);
    MB_birth(i).beta = t;
    MB_birth(i).epsilon = t;
end

