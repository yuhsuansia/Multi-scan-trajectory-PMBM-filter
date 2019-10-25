function [unknownPPP,trajectoryUpdMBM,trajectoryNewMBM] = ...
    updateStep(unknownPPP,trajectoryMBM,model,z,t)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENT

% Extract parameters from model
Pd = model.Pd;
H = model.H;
R = model.R;
lambda_fa = model.lambda_fa;

% Interpret sizes from inputs
% If it is a hypothesis with zero existence probability, no update
% e.g. n_invalid, n_valid, nupd = n_invalid + n_valid*(m+1);
nu = length(unknownPPP);
[~,m] = size(z);
n = length(trajectoryMBM);  % number of single target hypotheses to be updated
% number of single target hypotheses after update
if n > 0
    nupd = length(find([trajectoryMBM.r]~=0))*m + n;
else
    nupd = 0;
end

% Allocate memory for existing tracks (single target hypotheses)
trajectoryUpdMBM = repmat(struct,nupd,1);

% Keep record ancestor information, i.e., the index of the single target
% hypothesis being updated, need be called in N-scan pruning

% Update existing tracks
iupd = 0;   % initiate updating index
insideGating = true(nupd,1);
usedMeas = false(m,1);
for i = 1:n
    % first determine whether the single target hypothesis has valid
    % existence probability or not
    if trajectoryMBM(i).r == 0 % this case corresponds to non-existence single target hypothesis
        iupd = iupd+1;
        trajectoryUpdMBM(iupd).c = 0;
        trajectoryUpdMBM(iupd).r = 0;
        trajectoryUpdMBM(iupd).xpre = zeros(4,1);
        trajectoryUpdMBM(iupd).Ppre = zeros(4,4);
        trajectoryUpdMBM(iupd).x = zeros(4,1);
        trajectoryUpdMBM(iupd).P = zeros(4,4);
        trajectoryUpdMBM(iupd).l = [trajectoryMBM(i).l;[t,0]];
        trajectoryUpdMBM(iupd).a = trajectoryMBM(i).a;
        trajectoryUpdMBM(iupd).beta = 0;
        trajectoryUpdMBM(iupd).epsilon = 0;
    else
        % Create missed detection hypothesis
        iupd = iupd+1;
        temp = 1-trajectoryMBM(i).r+trajectoryMBM(i).r*(1-Pd*trajectoryMBM(i).w(end));
        trajectoryUpdMBM(iupd).c = trajectoryMBM(i).c-log(temp);
        trajectoryUpdMBM(iupd).w = [trajectoryMBM(i).w(1:end-1) trajectoryMBM(i).w(end)*(1-Pd)]/(1-trajectoryMBM(i).w(end)*Pd);
        trajectoryUpdMBM(iupd).r = trajectoryMBM(i).r*(1-Pd*trajectoryMBM(i).w(end))/temp;
        trajectoryUpdMBM(iupd).xpre = trajectoryMBM(i).xpre;
        trajectoryUpdMBM(iupd).Ppre = trajectoryMBM(i).Ppre;
        trajectoryUpdMBM(iupd).x = trajectoryMBM(i).x;
        trajectoryUpdMBM(iupd).P = trajectoryMBM(i).P;
        % If it is missed detection, add 0 to measurement history
        trajectoryUpdMBM(iupd).l = [trajectoryMBM(i).l;[t,0]];
        trajectoryUpdMBM(iupd).a = trajectoryMBM(i).a;
        trajectoryUpdMBM(iupd).beta = trajectoryMBM(i).beta;
        trajectoryUpdMBM(iupd).epsilon = trajectoryMBM(i).epsilon;
        
        % Create hypotheses with measurement updates
        S = H*trajectoryMBM(i).P(:,:,end)*H' + R;
        sqrt_det2piS = sqrt(det(2*pi*S));
        K = trajectoryMBM(i).P(:,:,end)*H'/S;
        Pplus = trajectoryMBM(i).P(:,:,end) - K*H*trajectoryMBM(i).P(:,:,end);
        
        for j = 1:m
            iupd = iupd+1;
            % do not update hypotheses with too small existence probability
            % at the current time step
            if trajectoryMBM(i).r*trajectoryMBM(i).w(end) < model.threshold
                insideGating(iupd) = false;
            else
                v = z(:,j) - H*trajectoryMBM(i).x(:,end);
                maha = v'/S*v;
                temp = exp(-0.5*maha)/sqrt_det2piS;
                % ellipsoidal gating
                if maha > model.gamma
                    insideGating(iupd) = false;
                else
                    usedMeas(j) = true;
                    trajectoryUpdMBM(iupd).c = trajectoryMBM(i).c-log(trajectoryMBM(i).r*Pd*trajectoryMBM(i).w(end)*temp);
                    trajectoryUpdMBM(iupd).w = 1;
                    trajectoryUpdMBM(iupd).r = 1;
                    trajectoryUpdMBM(iupd).xpre = trajectoryMBM(i).xpre;
                    trajectoryUpdMBM(iupd).Ppre = trajectoryMBM(i).Ppre;
                    trajectoryUpdMBM(iupd).x = trajectoryMBM(i).x;
                    trajectoryUpdMBM(iupd).x(:,end) = trajectoryMBM(i).x(:,end) + K*v;
                    trajectoryUpdMBM(iupd).P = trajectoryMBM(i).P;
                    trajectoryUpdMBM(iupd).P(:,:,end) = Pplus;
                    % Otherwise, add the index of the measurement at current scan
                    trajectoryUpdMBM(iupd).l = [trajectoryMBM(i).l;[t,j]];
                    trajectoryUpdMBM(iupd).a = trajectoryMBM(i).a;
                    trajectoryUpdMBM(iupd).beta = trajectoryMBM(i).beta;
                    trajectoryUpdMBM(iupd).epsilon = trajectoryMBM(i).epsilon(end);
                end
            end
        end
    end
end

unusedMeas = ~usedMeas;

% Only keep updated single target hypotheses
idx_keep = insideGating;
trajectoryUpdMBM = trajectoryUpdMBM(idx_keep);

% Allocate memory for new tracks, each new track contains two single target
% hypothese
trajectoryNewMBM = struct('w',0,'r',0,'xpre',zeros(4,1),'x',zeros(4,1),...
    'Ppre',zeros(4,4,1),'P',zeros(4,4,1),'l',[],'c',0,'a',[],'beta',0,'epsilon',0);
trajectoryNewMBM = repmat(trajectoryNewMBM,2*m,1);

% Allocate temporary memory for new tracks
Sk = zeros(2,2,nu);
Kk = zeros(4,2,nu);
Pk = zeros(4,4,nu);
ck = zeros(nu,1);
sqrt_det2piSk = zeros(nu,1);
yk = zeros(4,nu);

% Create a new track for each measurement by updating PPP with measurement
for k = 1:nu
    Pu = unknownPPP(k).P(:,:,end);
    Sk(:,:,k) = H*Pu*H' + R;
    sqrt_det2piSk(k) = sqrt(det(2*pi*Sk(:,:,k)));
    Kk(:,:,k) = Pu*H'/Sk(:,:,k);
    Pk(:,:,k) = Pu - Kk(:,:,k)*H*Pu;
end
for j = 1:m
    for k = 1:nu
        v = z(:,j) - H*unknownPPP(k).x(:,end);
        ck(k) = unknownPPP(k).lambda*Pd*exp(-0.5*v'/Sk(:,:,k)*v)/sqrt_det2piSk(k);
        yk(:,k) = unknownPPP(k).x(:,end) + Kk(:,:,k)*v;
    end
    C = sum(ck);
    % first single target hypothesis for measurement associated to previous
    % track, second for new track
    trajectoryNewMBM(2*j).w = 1;
    trajectoryNewMBM(2*j-1).c = 0;
    trajectoryNewMBM(2*j).c = -log(C + lambda_fa);
    trajectoryNewMBM(2*j).r = C/(C + lambda_fa);
    % only keep the most likely Gaussian component
    [~,maxidx] = max(ck);
    trajectoryNewMBM(2*j).xpre = unknownPPP(maxidx).x(:,end);
    trajectoryNewMBM(2*j).Ppre = unknownPPP(maxidx).P(:,:,end);
    ck = ck/C;
    trajectoryNewMBM(2*j).x = yk*ck;
    trajectoryNewMBM(2*j).P = zeros(4,4);
    for k = 1:nu
        v = trajectoryNewMBM(2*j).x - yk(:,k);
        trajectoryNewMBM(2*j).P = trajectoryNewMBM(2*j).P + ck(k)*(Pk(:,:,k) + v*v');
    end
    % for trajectory purpose
    trajectoryNewMBM(2*j-1).l = [t,0];    % add 0, if there is no new target
    trajectoryNewMBM(2*j).l = [t,j];      % otherwise, add measurement index
    % update track index
    if nupd==0
        trajectoryNewMBM(2*j-1).a = j;
        trajectoryNewMBM(2*j).a = j;
    else
        trajectoryNewMBM(2*j-1).a = trajectoryUpdMBM(end).a+j;
        trajectoryNewMBM(2*j).a = trajectoryUpdMBM(end).a+j;
    end
    trajectoryNewMBM(2*j-1).beta = 0;
    trajectoryNewMBM(2*j-1).epsilon = 0;
    trajectoryNewMBM(2*j).beta = t;
    trajectoryNewMBM(2*j).epsilon = t;
end

% If a newly created trajectory has too small existence probability,
% and the measurement used to create this trajectory does not fall inside
% the gates of any pre-existing tracks, remove it
idx_remain = true(2*m,1);
for i = 1:m
    if trajectoryNewMBM(2*i).r < model.threshold && unusedMeas(i)==true
        idx_remain(2*i-1:2*i) = false;
    end
end
trajectoryNewMBM = trajectoryNewMBM(idx_remain);

% Update (i.e., thin) intensity of unknown targets
for i = 1:nu
    unknownPPP(i).lambda = (1-Pd)*unknownPPP(i).lambda;
end

% Truncate low weight components
ss = [unknownPPP.lambda] > model.threshold;
unknownPPP = unknownPPP(ss);

end

