function [model,measlog,xlog,X] = gentruth(slideWindow,L,Pd,lfai,numtruth,Pmid,pmi)
%GENTRUTH: GENERATE SIMULATION TRUTH

if (pmi == 1) % in this case, targets are present to begin with
    birthtime = zeros(1,numtruth);
    deathtime = 81*ones(1,numtruth);
else % in other cases, a new target appears every 10 time steps
    birthtime = 10*(0:numtruth-1);
    deathtime = 51+10*(0:numtruth-1);
end

% Initialise model
T = 1;  % sampling time
% nearly constant velocity model
model.F = kron(eye(2),[1 T; 0 1]);
model.Q = 0.04*kron(eye(2),[T^3/3 T^2/2; T^2/2 T]);
model.Qc = chol(model.Q);
% linear measurement model
model.H = kron(eye(2),[1 0]);
model.R = 0.04*eye(2);
model.Rc = chol(model.R);
model.Pd = Pd;              % detection probability
model.Ps = 0.99;            % survival probability
model.existThresh = 0.4;    % existence probability used to extract estimates
model.threshold = 1e-4;     % pruning threshold     
model.slideWindow = slideWindow;      % N-scan pruning
model.L = L;                % L-scan trajectory density approximation

% Gating parameters
P_G= 0.999;                    %gate size in percentage
model.gamma= chi2inv(P_G,2);   %inv chi^2 dn gamma value

% Initialise new target parameter structure
model.lambdab = 0.05;           % expect one new target to arrive every 20 scans on average
model.xb = zeros(4,1);
model.Pb = diag([100 1 100 1].^2);

volume = 200*200;
model.lfai = lfai;              % expected number of false alarms (integral of lambda_fa)
model.lambda_fa = lfai/volume;  % intensity = expected number / state space volume

simlen = 81; % must be odd
midpoint = (simlen+1)/2;
numfb = midpoint-1;
xlog = cell(simlen,1);

X.tVec = birthtime + 1;
X.iVec = deathtime - birthtime;
X.xState = zeros(2,simlen,numtruth);

% Initialise at time midpoint and propagate forward and backwards

% Run forward and backward simulation process
for i = 1:numtruth
    x = chol(Pmid)'*randn(size(model.F,1),1);
    xf = x; xb = x;
    X.xState(:,midpoint,i) = x([1,3]);
    xlog{midpoint} = [xlog{midpoint} x];
    for t = 1:numfb
        xf = model.F*xf + model.Qc'*randn(size(model.F,1),1);
        xb = model.F\(xb + model.Qc'*randn(size(model.F,1),1));
        xlog{midpoint-t} = [xlog{midpoint-t} xb(:,midpoint-t>birthtime(i))];
        xlog{midpoint+t} = [xlog{midpoint+t} xf(:,midpoint+t<=deathtime(i))];
        if midpoint-t>birthtime(i)
            X.xState(:,midpoint-t,i) = xb([1,3],midpoint-t>birthtime(i));
        end
        if midpoint+t<=deathtime(i)
            X.xState(:,midpoint+t,i) = xf([1,3],midpoint+t<=deathtime(i));
        end
    end
end

measlog = cell(simlen,1);
for j = 1:simlen
    measlog{j} =  makemeas(xlog{j},model);
end

    function z = makemeas(x,model)
        % Generate target measurements (for every target)
        z = model.H*x + model.Rc'*randn(size(model.H,1),size(x,2));
        % Simulate missed detection process
        z = z(:,rand(size(z,2),1) < model.Pd);
        % Generate false alarms (spatially uniform on [-100,100]^2
        z = [z, 200*rand(size(model.H,1),poissrnd(model.lfai))-100];
        % Shuffle order
        z = z(:,randperm(size(z,2)));
    end

end


