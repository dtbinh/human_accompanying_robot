function outPara = pathPlannerGreedy(agent,inPara)
% define input arguments
x_h = inPara.pre_traj(:,2); % predicted human trajectory at one step ahead
safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;
h_v = inPara.h_v;
obs_info = inPara.obs_info;
safe_marg = inPara.safe_marg;

% define parameters
% non_intersect_flag = 0; % flag for showing whether imposing the non-intersection constraint
dt = 0.05; % time interval for sampling the points on the line of the robot's path. used for imposing non-intersection constriant
safe_marg2 = 0.1; % margin for the robot's path line from the obstacle

x_r = agent.currentPos(1:2);
r_v = agent.currentV;
maxA = agent.maxA;
% determine robot heading and next position
theta = calAngle(x_r-x_h);
x_r_next = x_r+r_v*[cos(theta);sin(theta)];
%{
% safety constraint
% determine if the waypoint is inside an obstacle
for jj = 1:size(obs_info,2)
%     if sum((x(1:2,ii+1)-obs_info(1:2,jj)).^2) - (obs_info(3,jj)+safe_marg)^2 < 0
        
    if non_intersect_flag == 1
        % line not intersecting with the obstacle
        n = floor(mpc_dt/dt);
        x0 = obs_info(1,jj); y0 = obs_info(2,jj);
        r = obs_info(3,jj);
        for kk = 0:n
            constr = [constr,sum((kk/n*x(1:2,ii+1)+(n-kk)/n*x(1:2,ii)-obs_info(1:2,jj)).^2)>=(r+safe_marg2)^2];
        end
    end
%     end
%}
% change robot speed to match the human's next speed
min_v = r_v - maxA*mpt_dt;
max_v = r_v + maxA*mpt_dt;
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

function angle = calAngle (vec)
len = norm(vec);
x = vec(1);
y = vec(2);
if x >= 0
    angle = asin(y/len);
else
    angle = pi - asin(y/len);
end
end