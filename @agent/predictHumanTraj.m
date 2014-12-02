function [outPara] = predictHumanTraj(agent,inPara)
% define input arguments
state = inPara.state;
hor = inPara.hor;
pre_type = inPara.pre_type;
dt = inPara.mpc_dt;

v = state([2,4]);
% initialization
pre_traj = zeros(2,hor+1); % current and predicted human future path
pre_traj(:,1) = state([1,3]);

if strcmp(pre_type,'extpol')
    for ii = 1:hor
        pre_traj(:,ii+1) = pre_traj(:,ii)+v*dt;
    end
elseif strcmp(pre_type,'IMM')
    
end

outPara = pre_traj;
end