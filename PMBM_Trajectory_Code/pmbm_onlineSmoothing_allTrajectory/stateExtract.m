function [est,multiTrajectories,Xest] = stateExtract(trajectoryEst)
% Extract target states

if isempty(trajectoryEst)
    est = zeros(4,0);
    multiTrajectories = struct('states',[],'beginTime',[],'endTime',[]);
    Xest.tVec = zeros(0,1);
    Xest.iVec = zeros(0,1);
    Xest.xState = zeros(2,0,0);
    return;
end

% MAP cardinality estimate
r = extractfield(trajectoryEst,'r');
% ensure that the case that r==1 counts
r(r==1) = 1-eps;                    
ss = false(size(r));
pcard = prod(1-r)*poly(-r./(1-r));
[~,n] = max(pcard);
[~,o] = sort(-r);
n = n - 1;
ss(o(1:n)) = true;

trajectoryEst = trajectoryEst(ss);
len = length(trajectoryEst);
est = zeros(4,len);

% extract the set of all trajectories

multiTrajectories = struct('states',[],'beginTime',[],'endTime',[]);
multiTrajectories = repmat(multiTrajectories,len,1);

idx = true(len,1);
for i = 1:len
    % select trajectory death time with the highest weight
    [~,tdeath] = max(trajectoryEst(i).w);
    % store estiamtes for set of current trajectories
    est(:,i) = trajectoryEst(i).x(:,end);
    if tdeath ~= length(trajectoryEst(i).w)
        idx(i) = false;
    end
    multiTrajectories(i).states = trajectoryEst(i).x(:,1:trajectoryEst(i).epsilon(tdeath)-trajectoryEst(i).beta+1);
    multiTrajectories(i).beginTime = trajectoryEst(i).beta;
    multiTrajectories(i).endTime = trajectoryEst(i).epsilon(tdeath);
end
% extract current set estimate
est = est(:,idx);

% re-construct output estimates for the use of trajectory metric
if isempty(multiTrajectories)
    Xest.tVec = zeros(0,1);
    Xest.iVec = zeros(0,1);
    Xest.xState = zeros(2,0,0);
else
    Xest.tVec = [multiTrajectories.beginTime];
    Xest.iVec = [multiTrajectories.endTime]-[multiTrajectories.beginTime]+1;
    Xest.xState = zeros(2,max([multiTrajectories.endTime]),length(Xest.tVec));
    for j = 1:length(Xest.tVec)
        Xest.xState(:,Xest.tVec(j):Xest.tVec(j)+Xest.iVec(j)-1,j) = multiTrajectories(j).states([1,3],:);
    end
end

end

