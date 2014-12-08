function outPara = pathPlanner(agent,inPara)
% include IPOPT in YALMIP
% addpath('D:\Program Files\MATLAB\2013a_crack\IPOPT3.11.8');
addpath('D:\Chang Liu\ipopt');
% define input arguments
x_h = inPara.pre_traj; % predicted human trajectory
hor = inPara.hor;
safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;
h_v = inPara.h_v;
obs_info = inPara.obs_info;
safe_marg = inPara.safe_marg;

% define parameters
non_intersect_flag = 0; % flag for showing whether imposing the non-intersection constraint
dt = 0.05; % time interval for sampling the points on the line of the robot's path. used for imposing non-intersection constriant
safe_marg2 = 0.1; % margin for the robot's path line from the obstacle
tmp_hor = hor;
r_hd = agent.currentPos(3);

while(tmp_hor > 0)
    % define MPC
    x = sdpvar(3,tmp_hor+1); %[x,y,v]
    u = sdpvar(2,tmp_hor); %[psi,a]
    
    % impose constraints
    % initial condition
    constr = [x(:,1) == [agent.currentPos(1:2);agent.currentV]];
    % constraints on future states
    inPara_cg = struct('hor',tmp_hor,'x',x,'u',u,'h_v',h_v,'mpc_dt',mpc_dt,...
        'safe_dis',safe_dis,'safe_marg',safe_marg,'x_h',x_h,'obs_info',obs_info,...
        'non_intersect_flag',non_intersect_flag,'obj',0,'constr',constr,...
        'agent',agent,'dt',dt,'safe_marg2',safe_marg2,'r_hd',r_hd);
    [obj,constr] = genMPC(inPara_cg); % generate constraints. contain a parameter that decides whether using the non-intersection constraints
    
    % solve MPC
    opt = sdpsettings('solver','ipopt','usex0',1,'debug',1);
    sol = optimize(constr,obj,opt);
    
    if sol.problem == 0
        opt_x = value(x); % current and future states
        opt_u = value(u); % future input
        if tmp_hor < hor
            opt_x = [opt_x,opt_x(:,end)*ones(1,hor-tmp_hor)]; % current and future states
            opt_u = [opt_u,zeros(size(opt_u,1),hor-tmp_hor)]; % future input
        end
        for ii = 1:tmp_hor
            for jj = 1:size(obs_info,2)
                % check if the line intersects with some obstacle
                n = floor(mpc_dt/dt);
%                 x0 = obs_info(1,jj); y0 = obs_info(2,jj);
                r = obs_info(3,jj);
                for kk = 0:n
                    tmp = sum((kk/n*opt_x(1:2,ii+1)+(n-kk)/n*opt_x(1:2,ii)-obs_info(1:2,jj)).^2) - (r+safe_marg2)^2;
                    if tmp < 0
                      non_intersect_flag = 1;
                      break
                    end
                end
                if tmp < 0
                    break
                end
            end
            if tmp < 0
                break
            end
        end
        if tmp >= 0
            break
        end
    else
        display('Fail to solve MPC')
        sol.info
        yalmiperror(sol.problem)
        tmp_hor = tmp_hor-1;
    end
end
outPara = struct('opt_x',opt_x,'opt_u',opt_u);
end

function dis = dis_pl(p,v1,v2) % point-to-line distance
% define line equation: ax+by+c = 0
a = v2(2)-v1(2);
b = v1(1)-v2(1);
c = v1(2)*v2(1)-v2(2)*v1(1);
% poin-to-line distance
dis = (a*p(1)+b*p(2)+c)^2/(a^2+b^2);
end

function [a,b,c] = getLine(v1,v2)
% define line equation: ax+by+c = 0
a = v2(2)-v1(2);
b = v1(1)-v2(1);
c = v1(2)*v2(1)-v2(2)*v1(1);
end

function [obj,constr] = genMPC(inPara)
hor = inPara.hor;
x = inPara.x;
u = inPara.u;
h_v = inPara.h_v;
mpc_dt = inPara.mpc_dt;
safe_dis = inPara.safe_dis;
safe_marg = inPara.safe_marg;
x_h = inPara.x_h;
obs_info = inPara.obs_info;
non_intersect_flag = inPara.non_intersect_flag;
obj = inPara.obj;
constr = inPara.constr;
agent = inPara.agent;
dt = inPara.dt;
safe_marg2 = inPara.safe_marg2;
r_hd = inPara.r_hd;

for ii = 1:hor
    hr_dis = sum((x(1:2,ii+1)-x_h(:,ii+1)).^2); % square of the Euclidean distance between human and robot
    if ii == 1
        obj = obj+hr_dis+0.1*(x(3,ii+1)-h_v)^2;%-0.05*log(hr_dis-safe_dis)+(sin(u(1,ii))-sin(r_hd))^2;%+0.5*u(2,ii)^2
    else
        obj = obj+hr_dis+0.1*(x(3,ii+1)-h_v)^2;%-0.05*log(hr_dis-safe_dis)+(sin(u(1,ii))-sin(u(1,ii-1)))^2;
    end
    % constraints on robot dynamics
    constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(3,ii)*[cos(u(1,ii));sin(u(1,ii))]*mpc_dt,...
        x(3,ii+1) == x(3,ii)+u(2,ii)*mpc_dt,x(3,ii+1)>=0,-agent.maxA<=u(2,ii)<=agent.maxA];   
    % constraint on safe distance
    constr = [constr,sum((x(1:2,ii+1)-x_h(:,ii+1)).^2) >= safe_dis^2];
    % constraint on obstacle avoidance
    % robot should not be inside the obstacles, i.e. robot waypoints should
    % not be inside the obstacle and the line connecting the waypoints 
    % should not intersect with the obstacle
%     [a,b,c] = getLine(x(1:2,ii+1),x(1:2,ii));
    for jj = 1:size(obs_info,2)
        % waypoints not inside the obstacle
        constr = [constr,sum((x(1:2,ii+1)-obs_info(1:2,jj)).^2) >= (obs_info(3,jj)+safe_marg)^2];
        if non_intersect_flag == 1
            % line not intersecting with the obstacle
            n = floor(mpc_dt/dt);
            x0 = obs_info(1,jj); y0 = obs_info(2,jj);
            r = obs_info(3,jj);
            for kk = 0:n
                constr = [constr,sum((kk/n*x(1:2,ii+1)+(n-kk)/n*x(1:2,ii)-obs_info(1:2,jj)).^2)>=(r+safe_marg2)^2];
            end
        end
    end    
end
end