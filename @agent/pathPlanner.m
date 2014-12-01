function outPara = pathPlanner(agent,inPara)
% define input arguments
x_h = inPara.pre_traj; % predicted human trajectory
hor = inPara.hor;
safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;

% define MPC
x = sdpvar(3,hor+1); %[x,y,v]
u = sdpvar(2,hor); %[psi,a]
obj = norm(x(1:2,end)-x_h(:,end),1);
constr = [];
for ii = 1:hor
    obj = obj+norm(x(1:2,ii)-x_h(:,ii),1);
    % constraints on robot dynamics
    constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(3,ii)*[cos(u(1,ii));sin(u(1,ii))]*mpc_dt,...
        x(3,ii+1) == x(3,ii)+u(2,ii)*mpc_dt,0<=u(2,ii)<=agent.maxA];   
    % constraint on safe distance
    constr = [constr,norm(x(1:2,ii+1)-x_h(:,ii+1),1) >= safe_dis];
    % constraint on obstacle avoidance
end
sol = optimize(constr,obj);
if sol.problem == 0
    fut_x = value(x); % current and future states
    fut_u = value(u); % future input
else
    display('Fail to solve MPC')
    sol.info
    yalmiperror(sol.problem)
end
outPara = struct('fut_x',fut_x,'fut_u',fut_u);