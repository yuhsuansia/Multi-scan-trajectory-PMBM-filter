function trajectoryMBM = approxSmoothedTrajectory(trajectoryMBM,trajectorySmoothEst,model)
% Use part of smoothed trajectories to replace unsmoothed ones

aMBM =[trajectoryMBM.a];
if ~isempty(trajectorySmoothEst)
    aSmoothEst = [trajectorySmoothEst.a];
else
    aSmoothEst = [];
end

numBernoulli = length(trajectoryMBM);
[lia,locb]  = ismember(aMBM,aSmoothEst);

for i = 1:numBernoulli
    if lia(i)==true && trajectoryMBM(i).r~=0
        trajectoryLength = size(trajectoryMBM(i).x,2);
        smoothingLength = min(trajectoryLength,model.slideWindow+model.L);
        if smoothingLength < trajectoryLength
            trajectoryMBM(i).x(:,end-smoothingLength+1) = trajectorySmoothEst(locb(i)).x(:,end-smoothingLength+1);
        end
    end
end

end