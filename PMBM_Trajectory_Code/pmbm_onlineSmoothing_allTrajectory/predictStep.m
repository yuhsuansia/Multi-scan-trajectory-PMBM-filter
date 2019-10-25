function [trajectoryMBM,unknownPPP] = predictStep(trajectoryMBM,unknownPPP,model,t)
%PREDICT: PREDICT MULTI-BERNOULLI AND POISSON COMPONENTS

% Get multi-Bernoulli prediction parameters from model
F = model.F;
Q = model.Q;
Ps = model.Ps;

% Interpret length of inputs
nb = length(model.lambdab);
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

nu = length(unknownPPP);
% Predict existing PPP intensity
for i = 1:nu
    unknownPPP(i).lambda = Ps*unknownPPP(i).lambda;
    unknownPPP(i).x = [unknownPPP(i).x F*unknownPPP(i).x(:,end)];
    unknownPPP(i).P = cat(3,unknownPPP(i).P,F*unknownPPP(i).P(:,:,end)*F'+Q);
    % beta remains the same, only update epsilon
    unknownPPP(i).epsilon = t;
end

% Incorporate birth intensity into PPP
for i = 1:nb
    unknownPPP(nu+i,1).lambda = model.lambdab(i);
    unknownPPP(nu+i,1).x = model.xb(:,i);
    unknownPPP(nu+i,1).P = model.Pb(:,:,i);
    unknownPPP(nu+i,1).beta = t;
    unknownPPP(nu+i,1).epsilon = t;
end



