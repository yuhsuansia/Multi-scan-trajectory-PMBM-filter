function [trajectoryEst,trajectoryUpdMBM] = dataAssoc(trajectoryUpdMBM,model)

% num of single target hypotheses
H = length(trajectoryUpdMBM);        

% n: num of tracks
aupd = [trajectoryUpdMBM.a];
unique_a = unique(aupd,'stable');
n = length(unique_a);

% number of single target hypotheses in track i, i belongs {1,...,n}
ns = zeros(n,1);
for i = 1:n
    ns(i) = length(find(aupd==unique_a(i)));
end

% calculate the length of trajectory of each single target hypothesis in
% tracks, used for determine the num of scans should be used
tralen = zeros(H,1);
for i = 1:H
    tralen(i) = size(trajectoryUpdMBM(i).l,1);
end

% cost of single target hypotheses
cupd = [trajectoryUpdMBM.c]';
c = cupd;

% construct binary indicator matrix for constraint (1): each track should
% only be used once; this constraint is neccessary for the implementation of dual
% decomposition, since sliding window is used.
A0 = zeros(n,H);
% for each track
idx = 0;
for i = 1:n
    A0(i,idx+1:idx+ns(i)) = 1;
    idx = idx+ns(i);
end

maxtralen = max(tralen); % maxtralen-scan data association
tratemp = cell(maxtralen,1);

% constriant equality
slideWindow = model.slideWindow; % apply sliding window
maxtralen(maxtralen>slideWindow) = slideWindow;

% construct binary indicator matrix for constraint (2): each measurement in
% each scan should only be used once
A = cell(maxtralen,1);
b = cell(maxtralen,1);

At = [];
bt = [];

for tl = 1:maxtralen
    % get number of single target hypotheses existed at current-tl+1 scan
    Htemp = length(find(tralen>=tl));
    % target trajectory of the current-tl+1 scan, tracks from left to
    % right, old to new
    tratemp{tl} = arrayfun(@(x) x.l(end-tl+1,2),trajectoryUpdMBM(tralen>=tl));
    % num of measurements in current-tl+1 scan, do not count missed detection and
    % non-exist new track
    measUnique = unique(tratemp{tl}(tratemp{tl}~=0),'stable');
    mtemp = length(measUnique);
    
    Atemp = false(mtemp,Htemp); % current-tl+1 scan
    for i = 1:mtemp
        Atemp(i,tratemp{tl}==measUnique(i)) = true;
    end
    Atemp = cat(2,Atemp,false(mtemp,H-Htemp));
    btemp = ones(mtemp,1);
    
    A{tl} = Atemp;
    b{tl} = btemp;
    
    At = cat(1,At,Atemp);
    bt = cat(1,bt,btemp);
end

Amatrix = [A0;At];
% branch and bound options
options = optimoptions('intlinprog','Display','off');

% pre-calculated parameters
measindices = cell(maxtralen,1);
idx_miss = cell(maxtralen,1);
is_preCalculated = cell(maxtralen,1);
tratl = cell(maxtralen,n);
miss = cell(maxtralen,n);

for tl = 1:maxtralen
    num_meas = size(A{tl},1);
    measindices{tl} = cell(num_meas,n);
    idx_miss{tl} = cell(num_meas,n);
    for j = 1:num_meas
        idx = 0;
        for i = 1:n
            % find single target hypotheses in track i that use this
            % measurement
            is_preCalculated{tl}{j,i} = find(A{tl}(j,idx+1:idx+ns(i))==true);
            if ~isempty(is_preCalculated{tl}{j,i})
                measindices{tl}{j,i} = tratemp{tl}(idx+1:idx+ns(i));
                idx_miss{tl}{j,i} = find(measindices{tl}{j,i}==0);
            end
            idx = idx+ns(i);
        end
    end
    
    idx = 0;
    for i = 1:n
        if  size(trajectoryUpdMBM(idx+1).l,1) >= tl
            tratl{tl}{i} = tratemp{tl}(idx+1:idx+ns(i));
            miss{tl}{i} = find(tratl{tl}{i}==0);
        end
        idx = idx+ns(i);
    end
end

%%
% if 1
% dual decomposition, solution is a binary indicator vector, decides which
% single target hypotheses are included in the "best" global hypotheses.

% subproblem t: min(c/t+\deltat)*u, s.t. [A0;At]*u = [b0;bt];

% Larange multiplier \delta is initialised with 0
delta = zeros(H,maxtralen);
% subproblem solutions
u_hat = false(H,maxtralen);

% initialise maximum num of iteration
numIteration = 0;
maxIteration = 1e1;
numRepetition = false;
% store the best feasible primal cost obtained so far (upper bound)
bestPrimalCost = inf;
uprimal = false(H,1);

while (numIteration<maxIteration&&~numRepetition)
    % get suboptimal solution for each subproblem
    subDualCost = zeros(maxtralen,1);
    for tl = 1:maxtralen
        % implementation using auction
        c_hat = c/maxtralen+delta(:,tl);
        % get number of measurements at scan tl
        num_meas = size(A{tl},1);
        % construct track to measurement assignment matrix at scan tl
        cost = ones(n,num_meas)*inf;
        % store assignment index of the single target hypothesis with the
        % minimum cost in each track
        idxCost = zeros(n,num_meas);
        for j = 1:num_meas
            idx = 0;
            % find single target hypotheses in track i that use this
            % measurement if found, find the single target hypothesis with 
            % the minimum cost, and record its index
            for i = 1:n
                if ~isempty(is_preCalculated{tl}{j,i})
                    [cost(i,j),idxmin] = min(c_hat(idx+is_preCalculated{tl}{j,i}));
                    if ~isempty(idx_miss{tl}{j,i})
                        cost(i,j) = cost(i,j) - min(c_hat(idx+idx_miss{tl}{j,i}));
                    end
                    idxCost(i,j) = idx+is_preCalculated{tl}{j,i}(idxmin);
                end
                idx = idx+ns(i);
            end
        end
        % find the most likely assignment using optimal assignment
        costInput = [cost inf*ones(n,n-num_meas)];
        costInput = costInput-min(costInput(:));
        [assignments,~] = assignmentoptimal(costInput);
        assignments(assignments>num_meas) = 0;
        indicator = false(H,1);
        for i = 1:n
            if assignments(i)>0
                indicator(idxCost(i,assignments(i))) = true;
            end
        end
        
        % if a track has no measurement assigned to it, choose the single
        % target hypotheses to be non-exist or miss if the track exists
        % before scan N-tl, (if exists after scan N-tl, of course, no 
        % measurment would be assigned)
        utemp = false(H,1);
        utemp(indicator) = true;
        
        idx = 0;
        for i = 1:n
            if ~any(utemp(idx+1:idx+ns(i)))
                if size(trajectoryUpdMBM(idx+1).l,1) >= tl+1
                    [~,idxtratl] = min(c_hat(miss{tl}{i}+idx));
                    utemp(idx+miss{tl}{i}(idxtratl)) = true;
                else
                    [~,idxmin] = min(c_hat(idx+1:idx+ns(i)));
                    utemp(idx+idxmin) = true;
                end
            end
            idx = idx+ns(i);
        end
        u_hat(:,tl) = utemp;
        subDualCost(tl) = c_hat'*u_hat(:,tl);
    end
    % change data from double to binary
    u_hat = logical(u_hat);
    u_hat_mean = sum(u_hat,2)/maxtralen;
    
    % All the subproblem solutions are equal means we have found the
    % optimal solution
    if isempty(find(u_hat_mean~=1&u_hat_mean~=0,1))
        uprimal = u_hat(:,1);
        break;
    end
    
    % calculate dual cost
    dualCosthat = sum(subDualCost);
    
    % find partial primal solution without conflicts
    idx_selectedHypo = u_hat_mean==1;
    idx_selectedHypo(tralen<=slideWindow) = false;
    idx_unselectedHypo = ~idx_selectedHypo;
    idx_uncertainMeasTracks = sum(Amatrix(:,idx_selectedHypo),2)==0;
    % if we are certain about the track, remove the single target
    % hypotheses of it
    for i = 1:n
        if idx_uncertainMeasTracks(i) == false
            idx_unselectedHypo(Amatrix(i,:)==1) = false;
        end
    end
    A_uncertain = Amatrix(idx_uncertainMeasTracks,idx_unselectedHypo);
    c_uncertain = c(idx_unselectedHypo);
    len_c_uncertain = length(c_uncertain);
    
    % implement branch and bound algorithm to find a feasible solution
    uprimal_uncertain = intlinprog(c_uncertain,1:len_c_uncertain,[],[],...
        sparse(A_uncertain),ones(size(A_uncertain,1),1),...
        zeros(len_c_uncertain,1),ones(len_c_uncertain,1),[],options);
    uprimal_uncertain = round(uprimal_uncertain);

    uprimalhat = u_hat_mean==1;
    uprimalhat(idx_unselectedHypo) = logical(uprimal_uncertain);
    
    % obtain primal cost
    bestPrimalCosthat = c'*uprimalhat;
    
    if bestPrimalCosthat < bestPrimalCost
        bestPrimalCost = bestPrimalCosthat;
        uprimal = uprimalhat;
        numRepetition = false;
    else
        % jump out the loop if the best primal cost obtained does not increase
        numRepetition = true;
    end
    
    % calculate relative gap between dual cost and primal cost
    gap = (bestPrimalCost - dualCosthat)/bestPrimalCost;
    if gap < 0.02
        % jump out the loop if the gap is too small
        break;
    end
    
    % calculate step size used in subgradient methods
    % calculate subgradient
    g = u_hat - u_hat_mean;
    % calculate step size used in subgradient method
    stepSize = (bestPrimalCost - dualCosthat)/(norm(g)^2);
    % update Lagrange multiplier
    delta = delta + stepSize*g;
    
    % increase index of iterations
    numIteration = numIteration+1;
end

u = uprimal;

%%

% single target hypotheses in the ML global association hypotheses updating
% pre-existing tracks
I = u(1:H)==1;
trajectoryEst = trajectoryUpdMBM(I);

% N-scan pruning
idx = 0;
idx_remain = false(H,1);
nc = slideWindow;
for i = 1:length(trajectoryEst)
    trajectoryLength = size(trajectoryEst(i).l,1);
    if trajectoryLength>=nc && ns(i)==1 && trajectoryEst(i).r==0
        % prune tracks with only null-hypothesis with length no less than nc
    else
        if trajectoryLength>=slideWindow
            traCompared = trajectoryEst(i).l(1:end-slideWindow+1,2);
            for j = idx+1:idx+ns(i)
                % check if measurement history information matches
                if isequal(trajectoryUpdMBM(j).l(1:end-slideWindow+1,2),traCompared)
                    idx_remain(j) = true;
                end
            end
        else
            idx_remain(idx+1:idx+ns(i)) = true;
        end
    end
    idx = idx+ns(i);
end

trajectoryUpdMBM = trajectoryUpdMBM(idx_remain);

% find single target hypotheses with existence probability smaller than a
% pre-defined threshold. Check if the pruning of such hypotheses will
% affact the solvability of the multi-scan assignment, if not, prune it
Aupd = At(:,1:H);
Aupd = Aupd(:,idx_remain);
rupd = [trajectoryUpdMBM.r]';
idx_smallExistenceProb = rupd<model.threshold;
idx_measSmallExistenceProb = sum(Aupd(:,idx_smallExistenceProb),2)>=1;
idx_highExistenceProb = rupd>=model.threshold;
idx_measHighExistenceProb = sum(Aupd(:,idx_highExistenceProb),2)>=1;
indicator = (idx_measSmallExistenceProb-idx_measHighExistenceProb)==1;
len = size(Aupd,2);
idx_remain = true(len,1);
idx_remain(sum(Aupd.*indicator)>0) = false;
trajectoryUpdMBM = trajectoryUpdMBM(idx_remain);

% remove tracks with only non-existence or low existence hypothesis
if ~isempty(trajectoryUpdMBM)
    aupd = [trajectoryUpdMBM.a]';
    Hupd = length(trajectoryUpdMBM);
    idx_remain = true(Hupd,1);
    for i = 1:Hupd
        if length(find(aupd==trajectoryUpdMBM(i).a))==1 && trajectoryUpdMBM(i).r < model.threshold
            idx_remain(i) = false;
        end
    end
    trajectoryUpdMBM = trajectoryUpdMBM(idx_remain);
end

end