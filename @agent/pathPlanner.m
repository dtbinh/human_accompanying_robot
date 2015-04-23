function outPara = pathPlanner(agent,inPara)
% include IPOPT in YALMIP
% addpath('D:\Program Files\MATLAB\2013a_crack\IPOPT3.11.8');
% addpath('D:\Chang Liu\ipopt');
addpath('/Users/changliu/Documents/MATLAB/Ipopt-3.11.8-linux64mac64win32win64-matlabmexfiles')
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
h_v_value = norm(h_v,2);

init_state = [agent.currentPos;agent.currentV];
while(tmp_hor > 0)
    % define MPC
    x = sdpvar(4,tmp_hor+1); %[x,y,theta,v]
    u = sdpvar(2,tmp_hor); %[w,a]
    
    % impose constraints
    % initial condition
    constr = [x(:,1) == init_state];
    % constraints on future states
    inPara_cg = struct('hor',tmp_hor,'x',x,'u',u,'h_v',h_v_value,'mpc_dt',mpc_dt,...
        'safe_dis',safe_dis,'safe_marg',safe_marg,'x_h',x_h,'obs_info',obs_info,...
        'non_intersect_flag',non_intersect_flag,'obj',0,'constr',constr,...
        'agent',agent,'dt',dt,'safe_marg2',safe_marg2,'init_state',init_state);
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
                % x0 = obs_info(1,jj); y0 = obs_info(2,jj);
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
        % safe_dis = safe_dis/2;   
        tmp_hor = tmp_hor-1;
    end
end

if tmp_hor == 0 % if the MPC fails, just find the input at the next step to maximize the humna-robot distance
    %{
    % define MPC
    x = sdpvar(4,hor+1); %[x,y,theta,v]
    u = sdpvar(2,hor); %[w,a]
    tmp_hor = 1;
    % impose constraints
    % initial condition
    constr = [x(:,1) == init_state];
    % constraints on future states
    %     for tmp_hor = 2:hor+1
    inPara_cg = struct('hor',tmp_hor,'x',x,'u',u,'h_v',h_v_value,'mpc_dt',mpc_dt,...
        'safe_dis',safe_dis,'safe_marg',safe_marg,'x_h',x_h,'obs_info',obs_info,...
        'non_intersect_flag',non_intersect_flag,'obj',0,'constr',constr,...
        'agent',agent,'dt',dt,'safe_marg2',safe_marg2,'init_state',init_state);
    [obj,constr] = genMPCinfea(inPara_cg); % generate constraints. contain a parameter that decides whether using the non-intersection constraints
    % solve MPC
    opt = sdpsettings('solver','ipopt','usex0',1,'debug',1);
    sol = optimize(constr,obj,opt);
    
    if sol.problem == 0
        opt_x = value(x); % current and future states
        opt_u = value(u); % future input
        opt_x(:,tmp_hor+2:end) = zeros(size(opt_x,1),size(opt_x,2)-tmp_hor-1);
        opt_u(:,tmp_hor+2:end) = zeros(size(opt_u,1),size(opt_u,2)-tmp_hor-1);
        %{
        for ii = 1:hor
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
        %}
    else
        display('Fail to solve MPC')
        sol.info
        yalmiperror(sol.problem)
        %         safe_dis = safe_dis/2;
        %         tmp_hor = tmp_hor-1;
    end
    %     end
    %}
    x_r = agent.currentPos(1:2);
    r_hd = agent.currentPos(3);
    r_v = agent.currentV;
    a_lb = agent.a_lb;
    a_ub = agent.a_ub;
    w_lb = agent.w_lb;
    w_ub = agent.w_ub;
    x_r_next = x_r+r_v*[cos(r_hd);sin(r_hd)]*mpc_dt;
    rh_dis_next = sqrt(sum((x_r_next - x_h(:,2)).^2));
    rh_dir = calAngle(x_h(:,2)-x_r_next); % direction from robot to human
    if (rh_dis_next >= safe_dis)
        % if robot will be outside of the collision region, then turn its
        % heading toward the human's next position
        min_hd = r_hd + w_lb*mpc_dt;
        max_hd = r_hd + w_ub*mpc_dt;
        if rh_dir<=max_hd && rh_dir>=min_hd
            r_hd_next = rh_dir;
        else
            hd_dif_min = min(abs(min_hd-rh_dir),abs(2*pi-abs(min_hd-rh_dir)));
            hd_dif_max = min(abs(max_hd-rh_dir),abs(2*pi-abs(max_hd-rh_dir)));
            if hd_dif_min < hd_dif_max
                r_hd_next = min_hd;
            else
                r_hd_next = max_hd;
            end
        end
        
        %     r_v_next = r_v + a_ub*mpc_dt;
    else
        % if robot will be inside collision region, then turn its
        % heading against the human's next position
        min_hd = r_hd + w_lb*mpc_dt;
        max_hd = r_hd + w_ub*mpc_dt;
        op_rh_dir = -rh_dir;
        op_rh_dir = op_rh_dir - floor(op_rh_dir/(2*pi))*2*pi; % opposite direction
        if op_rh_dir<=max_hd && op_rh_dir>=min_hd
            r_hd_next = op_rh_dir;
        else
            hd_dif_min = min(abs(min_hd-rh_dir),abs(2*pi-abs(min_hd-rh_dir)));
            hd_dif_max = min(abs(max_hd-rh_dir),abs(2*pi-abs(max_hd-rh_dir)));
            if hd_dif_min < hd_dif_max
                r_hd_next = max_hd;
            else
                r_hd_next = min_hd;
            end
            %     elseif op_rh_dir<min_hd
            %         r_hd_next = max_hd;
            %     elseif op_rh_dir>max_hd
            %         r_hd_next = min_hd;
        end
        
        %     r_v_next = r_v + a_lb*mpc_dt;
    end
    tmp = r_hd_next;
    tmp = tmp - 2*pi*floor(tmp/(2*pi));
    r_hd_next = tmp;


    % chang robot speed to match human's current estimated speed
    min_v = r_v + a_lb*mpc_dt;
    max_v = r_v + a_ub*mpc_dt;
    
    if (rh_dis_next >= 2*safe_dis)
        r_v_next = max_v;
    elseif (rh_dis_next >= safe_dis) && (rh_dis_next < 2*safe_dis)
        if norm(h_v,2) >= max_v
            r_v_next = max_v;
        elseif norm(h_v,2) <= min_v
            r_v_next = min_v;
        else
            r_v_next = norm(h_v,2);
        end
    else
        r_v_next = min_v;
    end
    r_v_next = max(r_v_next,0);
    
    opt_x = [[x_r;r_hd;r_v],[x_r_next*ones(1,hor);r_hd_next*ones(1,hor);r_v_next,zeros(1,hor-1)]];
    opt_u = [[(r_hd_next-r_hd)/mpc_dt;(r_v_next-r_v)/mpc_dt],zeros(2,hor-1)];
end

% normalize the heading to [0,2*pi)
for ii = 1:size(opt_x,2)
    tmp = opt_x(3,ii);
    tmp = tmp - floor(tmp/(2*pi))*2*pi;
    opt_x(3,ii) = tmp;
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
% init_state = inPara.init_state;

% [A,B,c] = linearize_model(init_state,mpc_dt);
for ii = 1:hor
    hr_dis = sum((x(1:2,ii+1)-x_h(:,ii+1)).^2); % square of the Euclidean distance between human and robot
    if ii == 1
        obj = obj+hr_dis+0.2*(x(4,ii+1)-h_v)^2;%-0.1*log(hr_dis-safe_dis^2);%+(sin(u(1,ii))-sin(r_hd))^2;%+0.5*u(2,ii)^2
    else
        obj = obj+hr_dis+0.2*(x(4,ii+1)-h_v)^2;%-0.1*log(hr_dis-safe_dis^2);%+(sin(u(1,ii))-sin(u(1,ii-1)))^2;
    end
    % constraints on robot dynamics
    constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(4,ii)*[cos(x(3,ii));sin(x(3,ii))]*mpc_dt,...
        x(3,ii+1) == x(3,ii) + u(1,ii)*mpc_dt, x(4,ii+1) == x(4,ii)+u(2,ii)*mpc_dt,...
       x(4,ii+1)>=0,agent.a_lb<=u(2,ii)<=agent.a_ub,agent.w_lb<=u(1,ii)<=agent.w_ub];%
%     constr = [constr,x(:,ii+1) == A*x(:,ii)+B*u(:,ii)+c];
%     constr = [constr,x(4,ii+1)>=0,agent.a_lb<=u(2,ii)<=agent.a_ub,-agent.maxW<=u(1,ii)<=agent.maxW];%
    % constraint on safe distance
    constr = [constr,sum((x(1:2,ii+1)-x_h(:,ii+1)).^2) >= safe_dis^2];
%     constr = [constr,max(x(1:2,ii+1)-x_h(:,ii+1)) >= safe_dis];
    % constraint on obstacle avoidance
    % robot should not be inside the obstacles, i.e. robot waypoints should
    % not be inside the obstacle and the line connecting the waypoints 
    % should not intersect with the obstacle
%     [a,b,c] = getLine(x(1:2,ii+1),x(1:2,ii));
%{
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
%}
end
end

function [obj,constr] = genMPCinfea(inPara)
% solve optimization problem when the initial condition is infeasible
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
init_state = inPara.init_state;

% [A,B,c] = linearize_model(init_state,mpc_dt);
% constr = [constr,sum((x(1:2,t)-x_h(:,t)).^2) >= safe_dis^2];
obj = -sum((x(1:2,2)-x_h(:,2)).^2);
for ii = 1:hor
    % constraints on robot dynamics
    constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(4,ii)*[cos(x(3,ii));sin(x(3,ii))]*mpc_dt,...
        x(3,ii+1) == x(3,ii) + u(1,ii)*mpc_dt, x(4,ii+1) == x(4,ii)+u(2,ii)*mpc_dt,...
       x(4,ii+1)>=0,agent.a_lb<=u(2,ii)<=agent.a_ub,agent.w_lb<=u(1,ii)<=agent.w_ub];
%     constr = [constr,x(:,ii+1) == A*x(:,ii)+B*u(:,ii)+c];
%     constr = [constr,x(4,ii+1)>=0,agent.a_lb<=u(2,ii)<=agent.a_ub,-agent.maxW<=u(1,ii)<=agent.maxW];%
    % constraint on safe distance
    
%     constr = [constr,max(x(1:2,ii+1)-x_h(:,ii+1)) >= safe_dis];
    % constraint on obstacle avoidance
    % robot should not be inside the obstacles, i.e. robot waypoints should
    % not be inside the obstacle and the line connecting the waypoints 
    % should not intersect with the obstacle
%     [a,b,c] = getLine(x(1:2,ii+1),x(1:2,ii));
%{
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
%}
end
end