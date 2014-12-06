function outPara = pathPlanner(agent,inPara)
% include IPOPT in YALMIP
addpath('D:\Program Files\MATLAB\2013a_crack\IPOPT3.11.8');
% define input arguments
x_h = inPara.pre_traj; % predicted human trajectory
hor = inPara.hor;
safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;
h_v = inPara.h_v;
obs_info = inPara.obs_info;
safe_marg = inPara.safe_marg;

while(hor > 0)
    % define MPC
    x = sdpvar(3,hor+1); %[x,y,v]
    u = sdpvar(2,hor); %[psi,a]
%     obj = 0;
    
    % impose constraints
    % initial condition
    constr = [x(:,1) == [agent.currentPos(1:2);agent.currentV]];
    % constraints on future states
    inPara_cg = struct('hor',hor,'x',x,'u',u,'h_v',h_v,'mpc_dt',mpc_dt,...
        'safe_dis',safe_dis,'safe_marg',safe_marg,'x_h',x_h,'obs_info',obs_info,...
        'non_intersect_flag',0,'obj',0,'constr',constr);
    [obj,constr] = genConstr(inPara_cg); % generate constraints. contain a parameter that decides whether using the non-intersection constraints
    
    % solve MPC
    opt = sdpsettings('solver','ipopt','usex0',1,'debug',1);
    sol = optimize(constr,obj,opt);
    
    % if MPC is solved, then check if non-intersection constraint is
    % violated. If so, add such constraint and redo the MPC. If can be
    % solved, then use this result; otherwise, decrease the horizon and
    % solve MPC again.
    if sol.problem == 0
        opt_x = value(x); % current and future states
        opt_u = value(u); % future input
        % check if violating the non-intersection constraint
        vio_flag = 0; % flag for violating the constraints. 1 means violation
        for ii = 1:hor
            [a,b,c] = getLine(opt_x(1:2,ii+1),opt_x(1:2,ii));
            for jj = 1:size(obs_info,2)
                x0 = obs_info(1,jj); y0 = obs_info(2,jj);
                r = obs_info(3,jj);
                tmp = (x0^2-r^2)*a^2+(y0^2-r^2)*b^2+c^2+2*a*b*x0*y0+2*a*c*x0+2*b*c*y0;
                if tmp < 0
                    vio_flag = 1; 
                    break
                end
            end
            if vio_flag == 1
                break
            end
        end
        if vio_flag == 1
            % initial condition
            constr = [x(:,1) == [agent.currentPos(1:2);agent.currentV]];
            % constraints on future states
            inPara_cg = struct('hor',hor,'x',x,'u',u,'h_v',h_v,'mpc_dt',mpc_dt,...
                'safe_dis',safe_dis,'safe_marg',safe_marg,'x_h',x_h,'obs_info',obs_info,...
                'non_intersect_flag',1,'obj',0,'constr',constr);
            [obj,constr] = genConstr(inPara_cg);
            % solve MPC
            opt = sdpsettings('solver','ipopt','usex0',1,'debug',1);
            sol = optimize(constr,obj,opt);
            if sol.problem == 0
                opt_x = value(x); % current and future states
                opt_u = value(u); % future input
            else
                display('Fail to solve MPC')
                sol.info
                yalmiperror(sol.problem)
                hor = hor - 1;
            end
        elseif vio_flag == 0
            break;
        end
    else
        display('Fail to solve MPC')
        sol.info
        yalmiperror(sol.problem)
        break
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

function [obj,constr] = genConstr(inPara)
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
for ii = 1:hor
    obj = obj+sum((x(1:2,ii+1)-x_h(:,ii+1)).^2)+0.1*u(2,ii)^2+0.1*(x(3,ii+1)-h_v)^2;
    % constraints on robot dynamics
    constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(3,ii)*[cos(u(1,ii));sin(u(1,ii))]*mpc_dt,...
        x(3,ii+1) == x(3,ii)+u(2,ii)*mpc_dt];%,-agent.maxA<=u(2,ii)<=agent.maxA];   
    % constraint on safe distance
    constr = [constr,sum((x(1:2,ii+1)-x_h(:,ii+1)).^2) >= safe_dis^2];
    % constraint on obstacle avoidance
    % robot should not be inside the obstacles, i.e. both robot waypoints
    % and the line connecting the waypoints should not intersect with the obstacle
    [a,b,c] = getLine(x(1:2,ii+1),x(1:2,ii));
    for jj = 1:size(obs_info,2)
        constr = [constr,sum((x(1:2,ii+1)-obs_info(1:2,jj)).^2) >= (obs_info(3,jj)+safe_marg)^2];
%         constr = [constr,dis_pl(obs_info(1:2,jj),x(1:2,ii+1),x(1:2,ii))>=obs_info(3,jj)^2];
        if non_intersect_flag == 1
            x0 = obs_info(1,jj); y0 = obs_info(2,jj);
            r = obs_info(3,jj);
            constr = [constr,(x0^2-r^2)*a^2+(y0^2-r^2)*b^2+c^2+2*a*b*x0*y0+2*a*c*x0+2*b*c*y0 >= 0.01];
        end
    end    
end
end