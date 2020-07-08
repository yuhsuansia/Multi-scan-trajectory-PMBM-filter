function trajectoryUpdMBM = updateStep(MB_birth,trajectoryMBM,model,z,t)
%UPDATE: CONSTRUCT COMPONENTS OF DISTRIBUTION UPDATED WITH MEASUREMENT

% Extract parameters from model
Pd = model.Pd;
H = model.H;
R = model.R;
lambda_fa = model.lambda_fa;

% Interpret sizes from inputs
% If it is a hypothesis with zero existence probability, no update
% e.g. n_invalid, n_valid, nupd = n_invalid + n_valid*(m+1);
nb = length(MB_birth);
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
% hypothesis being updated, can be used in N-scan pruning

% Update existing tracks
iupd = 0;   % initiate updating index
insideGating = true(nupd,1);
usedMeas = false(m,1);
for i = 1:n
    % first determine whether the single target hypothesis has valid
    % existence probability or not
    if trajectoryMBM(i).r == 0 % this case corresponds to non-existence single target hypothesis
        iupd = iupd+1;
        trajectoryUpdMBM(iupd).c = trajectoryMBM(i).c;
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
        trajectoryUpdMBM(iupd).c = trajectoryMBM(i).c-log(1-Pd*trajectoryMBM(i).w(end));
        trajectoryUpdMBM(iupd).w = [trajectoryMBM(i).w(1:end-1) trajectoryMBM(i).w(end)*(1-Pd)]/(1-trajectoryMBM(i).w(end)*Pd);
        trajectoryUpdMBM(iupd).r = 1;
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
                if maha > model.gamma % ellipsoidal gating
                    insideGating(iupd) = false;
                else
                    usedMeas(j) = true;
                    trajectoryUpdMBM(iupd).c = trajectoryMBM(i).c-log(Pd*trajectoryMBM(i).w(end)*temp);
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

% Create a new track for each Bernoulli component in MB birth, each new
% track contains (#measurements+2) single trajectory hypotheses, one for
% not born, one for missed detection, #measurements for measurement updates
trajectoryNewMBM = struct('w',1,'r',0,'xpre',zeros(4,1),'x',zeros(4,1),...
    'Ppre',zeros(4,4,1),'P',zeros(4,4,1),'l',[],'c',0,'a',[],'beta',0,'epsilon',0);
trajectoryNewMBM = repmat(trajectoryNewMBM,(m+2)*nb,1);

% Allocate temporary working for new tracks
Sk = zeros(2,2,nb);
Kk = zeros(4,2,nb);
Pk = zeros(4,4,nb);
sqrt_det2piSk = zeros(nb,1);

idx_remain = true((m+2)*nb,1);

for i = 1:nb
    Pb = MB_birth(i).P;
    Sk(:,:,i) = H*Pb*H' + R;
    sqrt_det2piSk(i) = sqrt(det(2*pi*Sk(:,:,i)));
    Kk(:,:,i) = Pb*H'/Sk(:,:,i);
    Pk(:,:,i) = Pb - Kk(:,:,i)*H*Pb;
    
    % not born hypothesis
    trajectoryNewMBM((i-1)*(m+2)+1).c = -log(1-MB_birth(i).r);
    trajectoryNewMBM((i-1)*(m+2)+1).l = [t,0];
    if nupd==0
        trajectoryNewMBM((i-1)*(m+2)+1).a = i;
    else
        trajectoryNewMBM((i-1)*(m+2)+1).a = trajectoryUpdMBM(end).a+i;
    end
    
    % missed detection hypothesis
    trajectoryNewMBM((i-1)*(m+2)+2).c = -log((1-Pd)*MB_birth(i).r);
    trajectoryNewMBM((i-1)*(m+2)+2).w = 1-Pd;
    trajectoryNewMBM((i-1)*(m+2)+2).r = 1;
    trajectoryNewMBM((i-1)*(m+2)+2).xpre = MB_birth(i).x;
    trajectoryNewMBM((i-1)*(m+2)+2).Ppre = MB_birth(i).P;
    trajectoryNewMBM((i-1)*(m+2)+2).x = MB_birth(i).x;
    trajectoryNewMBM((i-1)*(m+2)+2).P = MB_birth(i).P;
    trajectoryNewMBM((i-1)*(m+2)+2).l = [t,0];
    trajectoryNewMBM((i-1)*(m+2)+2).a = trajectoryNewMBM((i-1)*(m+2)+1).a;
    trajectoryNewMBM((i-1)*(m+2)+2).beta = MB_birth(i).beta;
    trajectoryNewMBM((i-1)*(m+2)+2).epsilon = MB_birth(i).epsilon;
    
    % measurement update
    for j = 1:m
        v = z(:,j) - H*MB_birth(i).x;
        ck = MB_birth(i).r*Pd*exp(-0.5*v'/Sk(:,:,i)*v)/sqrt_det2piSk(i);
        % If a newly created trajectory with has too small existence probability,
        % and the measurement used to create this trajectory does not fall inside
        % the gates of any pre-existing tracks, remove it
        if ck/(ck+lambda_fa) < model.threshold && unusedMeas(j)==true
            idx_remain((i-1)*(m+2)+2+j) = false;
        else
            usedMeas(j) = true;
            trajectoryNewMBM((i-1)*(m+2)+2+j).c =-log(ck);
            trajectoryNewMBM((i-1)*(m+2)+2+j).r = 1;
            trajectoryNewMBM((i-1)*(m+2)+2+j).xpre = MB_birth(i).x;
            trajectoryNewMBM((i-1)*(m+2)+2+j).Ppre = MB_birth(i).P;
            trajectoryNewMBM((i-1)*(m+2)+2+j).x = MB_birth(i).x + Kk(:,:,i)*v;
            trajectoryNewMBM((i-1)*(m+2)+2+j).P = Pk(:,:,i);
            % Otherwise, add the index of the measurement at current scan
            trajectoryNewMBM((i-1)*(m+2)+2+j).l = [t,j];
            trajectoryNewMBM((i-1)*(m+2)+2+j).a = trajectoryNewMBM((i-1)*(m+1)+2).a;
            trajectoryNewMBM((i-1)*(m+2)+2+j).beta = MB_birth(i).beta;
            trajectoryNewMBM((i-1)*(m+2)+2+j).epsilon = t;
        end
    end
end

trajectoryNewMBM = trajectoryNewMBM(idx_remain);

% Create a track for each used measurement, each track contains two single
% trajectory hypothesis, one representing clutter, one representing the
% measurement is associated to one of the pre-existing tracks
trajectoryClutterMBM = struct('w',0,'r',0,'xpre',zeros(4,1),'x',zeros(4,1),...
    'Ppre',zeros(4,4,1),'P',zeros(4,4,1),'l',[],'c',0,'a',[],'beta',0,'epsilon',0);
trajectoryClutterMBM = repmat(trajectoryClutterMBM,2*m,1);
for j = 1:m
    trajectoryClutterMBM(2*j-1).l = [t,0];
    trajectoryClutterMBM(2*j).l = [t,j];
    trajectoryClutterMBM(2*j-1).a = trajectoryNewMBM(end).a+j;
    trajectoryClutterMBM(2*j).a = trajectoryNewMBM(end).a+j;
    trajectoryClutterMBM(2*j-1).c = 0;
    trajectoryClutterMBM(2*j).c = -log(lambda_fa);
end
trajectoryClutterMBM = trajectoryClutterMBM(logical(kron(usedMeas,[1;1])));

if nupd==0
    trajectoryUpdMBM = [trajectoryNewMBM;trajectoryClutterMBM];
else
    trajectoryUpdMBM = [trajectoryUpdMBM;trajectoryNewMBM;trajectoryClutterMBM];
end

