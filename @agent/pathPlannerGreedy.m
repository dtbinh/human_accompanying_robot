function outPara = pathPlannerGreedy(agent,inPara)
% define input arguments
x_h = inPara.pre_traj; % predicted human trajectory
hor = inPara.hor;
safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;
h_v = inPara.h_v;
obs_info = inPara.obs_info;
safe_marg = inPara.safe_marg;
guess_x = inPara.guess_x;
guess_u = inPara.guess_u;
cmft_dis = inPara.cmft_dis;

% define parameters
non_intersect_flag = 0; % flag for showing whether imposing the non-intersection constraint
dt = 0.05; % time interval for sampling the points on the line of the robot's path. used for imposing non-intersection constriant
safe_marg2 = 0.1; % margin for the robot's path line from the obstacle
h_v_value = norm(h_v,2);

init_state = [agent.currentPos;agent.currentV];

% define MPC
x = sdpvar(4,hor+1); %[x,y,theta,v]
u = sdpvar(2,hor); %[w,a]

% impose constraints
% initial condition
constr = [x(:,1) == init_state];
% constraints on future states
inPara_cg = struct('hor',hor,'x',x,'u',u,'h_v',h_v_value,'mpc_dt',mpc_dt,...
    'safe_dis',safe_dis,'safe_marg',safe_marg,'x_h',x_h,'obs_info',obs_info,...
    'non_intersect_flag',non_intersect_flag,'obj',0,'constr',constr,...
    'agent',agent,'dt',dt,'safe_marg2',safe_marg2,'init_state',init_state,...
    'cmft_dis',cmft_dis);
[obj,constr] = genMPC(inPara_cg); % generate constraints. contain a parameter that decides whether using the non-intersection constraints

% in this mpc part, the strategy for dealing with infeasibility come from
% the pathPlanner.m in human_robot_search
while(1)
    % solve MPC
    if ~isempty(guess_x)
        assign(x,guess_x); % assign initial solution
        assign(u,guess_u); % assign initial input
    end
    
    opt = sdpsettings('solver','snopt','usex0',1,'debug',1);
    sol = optimize(constr,obj,opt);
    
    if sol.problem == 0
        disp('NLP solved successfully')
        opt_x = value(x); % current and future states
        opt_u = value(u); % future input
        break
        % temporarily remove the collision check. will reactivate this part
        % later
        %{
        for ii = 1:hor
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
        %}
    else
        display('Fail to solve MPC, will change solver')
        sol.info
        yalmiperror(sol.problem)
        % safe_dis = safe_dis/2;   
        figure % use empty figure to show that this if statement is executed
        
        % from the pathPlanner.m in human_robot_search:
        % follow the MATLAB's suggestion: find a feasible point and use as
        % the initial solution for the original problem
        tmp_obj1 = obj;
        assign(x,guess_x); % assign initial solution
        assign(u,guess_u); % assign initial input
%         tmp_opt1 = sdpsettings('solver','snopt','usex0',1,'debug',1,'verbose',1);
        tmp_opt1 = sdpsettings('solver','fmincon','usex0',1,'debug',1,'verbose',1);
        tmp_opt1.Algorithm = 'interior-point';
        sol = optimize(constr,tmp_obj1,tmp_opt1);
        if sol.problem ~= 0
            % if cannot find the feasible point using fmincon, try 0 as the
            % objective function
            display('Fail to solve the original MPC using NLP, will use NLP with 0 obj')  
            tmp_obj2 = 0;
            assign(x,guess_x); % assign initial solution
            assign(u,guess_u); % assign initial input
            tmp_opt2 = sdpsettings('solver','fmincon','usex0',1,'debug',1,'verbose',1);
            tmp_opt2.Algorithm = 'sqp';
            sol = optimize(constr,tmp_obj2,tmp_opt2);
            
            if sol.problem ~= 0
                display('Fail to solve the original MPC using NLP, will use reactive method')
                %{
                display('Fail to solve the original MPC using NLP, will use LP with 0 obj')
                % linearize system
                [A,B,c] = agent.linearize_model([agent.currentPos(1:3);agent.currentV],mpc_dt);
                tmp_constr3 = [x(:,1) == init_state];
                for ii = 1:hor
                    tmp_constr3 = [tmp_constr3,x(:,ii+1) == A*x(:,ii)+B*u(:,ii)+c,...
                        agent.w_lb<=u(1,ii)<=agent.w_ub,agent.a_lb<=u(2,ii)<=agent.a_ub,...
                        x(1:2,ii+1) >= 0];
                end
                tmp_obj3 = 0;
                tmp_opt3 = sdpsettings('solver','linprog','debug',1,'verbose',0);
                sol = optimize(tmp_constr3,tmp_obj3,tmp_opt3);
                
                if sol.problem ~= 0
                %}
                    % if the MPC fails, just find the input at the next 
                    % step to maximize the humna-robot distance
                    x_r = agent.currentPos(1:2);
                    r_hd = agent.currentPos(3);
                    r_v = agent.currentV;
                    a_lb = agent.a_lb;
                    a_ub = agent.a_ub;
                    w_lb = agent.w_lb;
                    w_ub = agent.w_ub;
                    x_r_next = x_r+r_v*[cos(r_hd);sin(r_hd)]*mpc_dt;
                    rh_dis_next = sqrt(sum((x_r_next - x_h).^2));
                    rh_dir = calAngle(x_h-x_r_next); % direction from robot to human
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
                        end
                    end
                    tmp = r_hd_next;
                    tmp = tmp - 2*pi*floor(tmp/(2*pi));
                    r_hd_next = tmp;
                    
                    % change robot speed to match human's current estimated speed
                    min_v = r_v + a_lb*mpc_dt;
                    max_v = r_v + a_ub*mpc_dt;
                    
                    if (rh_dis_next >= 1.2*safe_dis)
                        r_v_next = max_v;
                    elseif ((rh_dis_next >= 0.8*safe_dis) && (rh_dis_next < 1.2*safe_dis))
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
                    guess_x = opt_x
                    guess_u = opt_u
                    break
                %{
                elseif sol.problem == 0
                    disp('LP with 0 obj solved successfully')
                    guess_x = value(x);
                    guess_u = value(u);
                    opt_x = guess_x
                    opt_u = guess_u
                    break
                end
                %}
            elseif sol.problem == 0
                disp('NLP with 0 obj solved successfully')
                guess_x = value(x)
                guess_u = value(u)
            end
        elseif sol.problem == 0
            disp('NLP solved successfully using another solver')
            guess_x = value(x)
            guess_u = value(u)
        end
    end
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
cmft_dis = inPara.cmft_dis;
% init_state = inPara.init_state;

% [A,B,c] = linearize_model(init_state,mpc_dt);
for ii = 1:hor
    % objective
    hr_dis = abs(sum((x(1:2,ii+1)-x_h).^2)-cmft_dis^2); % square of the Euclidean distance between human and robot
    if ii == 1
        obj = obj+hr_dis+2*(x(4,ii+1)-h_v)^2;%-0.1*log(hr_dis-safe_dis^2);%+(sin(u(1,ii))-sin(r_hd))^2;%+0.5*u(2,ii)^2
    else
        obj = obj+hr_dis+2*(x(4,ii+1)-h_v)^2;%-0.1*log(hr_dis-safe_dis^2);%+(sin(u(1,ii))-sin(u(1,ii-1)))^2;
    end
    
    % constraints
    % constraints on robot dynamics
    constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(4,ii)*[cos(x(3,ii));sin(x(3,ii))]*mpc_dt,...
        x(3,ii+1) == x(3,ii) + u(1,ii)*mpc_dt, x(4,ii+1) == x(4,ii)+u(2,ii)*mpc_dt,...
       x(4,ii+1)>=0,agent.a_lb<=u(2,ii)<=agent.a_ub,agent.w_lb<=u(1,ii)<=agent.w_ub];%
%     constr = [constr,x(:,ii+1) == A*x(:,ii)+B*u(:,ii)+c];
%     constr = [constr,x(4,ii+1)>=0,agent.a_lb<=u(2,ii)<=agent.a_ub,-agent.maxW<=u(1,ii)<=agent.maxW];%
    
    % constraint on safe distance between human and robot
    constr = [constr,sum((x(1:2,ii+1)-x_h).^2) >= safe_dis^2];
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