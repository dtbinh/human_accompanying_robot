function [outPara] = predictHumanTraj(agent,inPara)
% define input arguments
state = inPara.state;
hor = inPara.hor;
pre_type = inPara.pre_type;

v = state([2,4]);
% initialization
pre_traj = zeros(2,hor+1); % predicted human future path, starting with his current position
pre_traj(:,1) = state([1,3]);

if strcmp(pre_type,'extpol')
    for ii = 1:hor
        pre_traj(:,ii+1) = pre_traj(:,ii)+v;
    end
elseif strcmp(pre_type,'IMM')
    
end

outPara = pre_traj;
end