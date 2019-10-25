function trajectorySmoothEst = backforwardFiltering(trajectoryEst,trajectorySmoothEstPre,model)
% Backforward filtering for trajectories in the most likely global
% hypothesis

% Input also include smoothing result of last time scan, try to reuse the
% calculated smoothing gain for the most likely global hypothesis

% Select trajectory with existence probability larger than 0
rest = [trajectoryEst.r];
trajectoryEst = trajectoryEst(rest>0);

% Preparation
trajectorySmoothEst = trajectoryEst;
numTrajectory = length(trajectoryEst);

if ~isempty(trajectorySmoothEstPre)
    apre = [trajectorySmoothEstPre.a];
    a = [trajectorySmoothEst.a];
    [lia,locb]  = ismember(a,apre);
    for i = 1:numTrajectory
        trajectoryLength = size(trajectoryEst(i).x,2);
        smoothingLength = min(trajectoryLength,model.slideWindow+model.L);
        trajectorySmoothEst(i).G = zeros(4,4,smoothingLength-1);
        if lia(i)==true && ~isempty(trajectorySmoothEstPre(locb(i)).G)
            l = trajectorySmoothEst(i).l(:,2);
            lpre = [trajectorySmoothEstPre(locb(i)).l(:,2);-1];
            idx = find((l-lpre)~=0,1);
            if idx > 1
                idxDiverge = min(smoothingLength,trajectoryLength-(idx-1))-1;
                trajectorySmoothEst(i).G(:,:,1:idxDiverge) = trajectorySmoothEstPre(locb(i)).G(:,:,1:idxDiverge);
            end
        end
    end
end

% Backward filtering
for i = 1:numTrajectory
    trajectoryLength = size(trajectoryEst(i).x,2);
    smoothingLength = min(trajectoryLength,model.slideWindow+model.L);
    
    if isempty(trajectorySmoothEstPre)
        trajectorySmoothEst(i).G = zeros(4,4,smoothingLength-1);
    end
    for j = 2:smoothingLength
        if ~any(trajectorySmoothEst(i).G(:,:,smoothingLength-j+1))
            trajectorySmoothEst(i).G(:,:,smoothingLength-j+1) = ...
                trajectoryEst(i).P(:,:,end-j+1)*model.F'/trajectoryEst(i).Ppre(:,:,end-j+2);
        end
        G = trajectorySmoothEst(i).G(:,:,smoothingLength-j+1);
        trajectorySmoothEst(i).x(:,end-j+1) = trajectoryEst(i).x(:,end-j+1) + ...
            G*(trajectorySmoothEst(i).x(:,end-j+2) - trajectoryEst(i).xpre(:,end-j+2));
        trajectorySmoothEst(i).P(:,:,end-j+1) = trajectoryEst(i).P(:,:,end-j+1) - ...
            G*(trajectoryEst(i).Ppre(:,:,end-j+2) - trajectorySmoothEst(i).P(:,:,end-j+2))*G';
    end
end

trajectorySmoothEst = rmfield(trajectorySmoothEst,{'xpre','Ppre'});

end