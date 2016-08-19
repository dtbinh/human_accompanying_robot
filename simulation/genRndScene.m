function [traj,traj_dense,traj_mpc,traj_mpc_dense] = genRndScene()
%%% randomly generate scenarios for simulation
% Chang Liu 8/6/16
% parameters to tune:
% s_pos: starting position of the human
% shape_set: the set of obstacles shapes
% obs_num: total number of obstacles
% tar_num: total number of targets
% max_spd_set: the maximum speed of the human between a pair of obstacles

clc;
close all
clear all

% basic setup
fld_cor = [0,0,100,100]; % size of the field [xmin,ymin,xmax,ymax]
s_pos = [5;5]; % starting position of the human

%% define obstacle
% define basic obstacle shapes
% {'r',[w,h]}: rectangle with horizontal size a (width) and vertical size b
% (height).
% {'c',a}: circle with radius a
shape_set = {{'r',[5,20]},{'r',[20,5]},{'r',[10,10]},{'c',5}}; 

% generate obstacle positions
obs_num = 5; % specifiy total number of obstacles
obs_set = cell(obs_num,1); % save info of obstacles
count = 0; % counter for obstacle

while (count < obs_num)
    % randomly generate obstacles until we have the numble of obstacle is obs_num
    shape_idx = randi(length(shape_set),1); % randomly choose a shape
    tmp_shape = shape_set{shape_idx};
    
    % generate position of an obstacle
    tmp_pos = [randi([fld_cor(1),fld_cor(3)],1);randi([fld_cor(2),fld_cor(4)],1)]; 
    if strcmp(tmp_shape{1},'c')        
    % for a circle, the generated position is its center
        radius = tmp_shape{2};
        
        % keep this obstacle if it falls within the field and do not
        % collide with other obstacles. otherwise, generate a new one
        if (tmp_pos(1)+radius <= fld_cor(3)) && (tmp_pos(1)-radius >= fld_cor(1)) ...
                && (tmp_pos(2)+radius <= fld_cor(4)) && (tmp_pos(2)-radius >= fld_cor(2))            
            % if inside the field, check collision then 
            new_shape = {shape_set{shape_idx}{1},shape_set{shape_idx}{2},tmp_pos};
            is_colliding = false;            
            for jj = 1:count                
                other_shape = obs_set{jj};                
                is_colliding = collisionCheck(new_shape,other_shape);
                if is_colliding
                    break;
                end
            end
            if is_colliding
                continue;
            end
            count = count + 1;
            obs_set{count} = new_shape;            
        end
    elseif strcmp(tmp_shape{1},'r')
        % for rectangle, the generated position is its lower-left coordinate
        w = tmp_shape{2}(1);
        h = tmp_shape{2}(2);
        
        % keep this obstacle if it falls within the field and do not
        % collide with other obstacles. otherwise, generate a new one
        if (tmp_pos(1)+w <= fld_cor(3)) && (tmp_pos(2)+h <= fld_cor(4))
            % if inside the field, check collision then   
            new_shape = {shape_set{shape_idx}{1},shape_set{shape_idx}{2},tmp_pos};
            is_colliding = false;
            for jj = 1:count
                other_shape = obs_set{jj};                
                is_colliding = collisionCheck(new_shape,other_shape);
                if is_colliding
                    break;
                end
            end
            if is_colliding
                continue
            end
            count = count + 1;
            obs_set{count} = new_shape;            
        end
    end
end

%% define targets locations
d_margin = 10; % we design the targets to lie within a smaller rectangular area than the whole field: [xmin+d_margin,ymin+d_margin,xmax-d_margin,ymax-d_margin]
tar_area = fld_cor+[d_margin,d_margin,-d_margin,-d_margin];

obs_margin = 2; % the target should be away from the boundaries of obstacles by this margin

% specify target number
tar_num = 5; 
tar_cnt = 1; % counter for target
tar_pos = zeros(2,tar_num); % save targets' positions

while (tar_cnt <= tar_num)
    re_gen = false; % if true, the generated point is inside an obstacle and new regeneration is needed
    
    % generate target position
    tmp_pos = [tar_area(3)-tar_area(1);tar_area(4)-tar_area(2)].*...
        [rand(1);rand(1)]+[tar_area(1);tar_area(2)];
    
    % check if the generated point is in an obstacle
    for ii = 1:obs_num 
        tmp_shape = obs_set{ii};
        if strcmp(tmp_shape{1},'c')
            % for circle, the generated position is its center
            radius = tmp_shape{2};
            if norm(tmp_pos-tmp_shape{3}) <= radius+obs_margin
                re_gen = true;
                break
            end
        elseif strcmp(tmp_shape{1},'r')
            % for rectangle, the generated position is its lower-left coordinate
%             w = tmp_shape{2}(1);
%             h = tmp_shape{2}(2);
%             ll_cor = tmp_shape{3}; % lower-left coordinate
            w = sqrt(2)*tmp_shape{2}(1);
            h = sqrt(2)*tmp_shape{2}(2);
            ll_cor = tmp_shape{3}-1/2*(sqrt(2)-1)*[w;h];
            
            if (tmp_pos(1) >= ll_cor(1)-obs_margin) && (tmp_pos(2) >= ll_cor(2)-obs_margin) && ...
                    (tmp_pos(1) <= ll_cor(1)+w+obs_margin) && (tmp_pos(2) <= ll_cor(2)+h+obs_margin)
                re_gen = true;
                break
            end
        end
    end
    
    if ~re_gen
        tar_pos(:,tar_cnt) = tmp_pos;
        tar_cnt = tar_cnt + 1;
    end
end
%}

%% draw the scenario
figure
hold on
rectangle('Position',[fld_cor(1),fld_cor(2),fld_cor(3)-fld_cor(1),fld_cor(4)-fld_cor(2)]);

% draw obstacles
delta_ang = linspace(0,2*pi,20);
for ii = 1:obs_num % check if the generated point is in an obstacle
    tmp_shape = obs_set{ii};
    if strcmp(tmp_shape{1},'c')
        % for circle, the generated position is its center
        radius = tmp_shape{2};
        center = tmp_shape{3};
        plot(center(1)+radius*cos(delta_ang),center(2)+radius*sin(delta_ang),'LineWidth',2);
        
    elseif strcmp(tmp_shape{1},'r')
        % for rectangle, the generated position is its lower-left coordinate
        w = tmp_shape{2}(1);
        h = tmp_shape{2}(2);
        ll_cor = tmp_shape{3}; % lower-left coordinate
        all_cor = [ll_cor,ll_cor+[w;0],ll_cor+[w;h],ll_cor+[0;h]];
        plot([all_cor(1,:),all_cor(1,1)],[all_cor(2,:),all_cor(2,1)],'LineWidth',2);
    end
end

% draw targets
for ii = 1:size(tar_pos,2)
    plot(tar_pos(1,ii),tar_pos(2,ii),'r*','MarkerSize',10)
end
%}

%% find feasible path between targets
%
% this code usese RRT* to find a path
hold on
addpath('rrtstar')

% spd_intv_set = [2.5 3;1.5 2;0.5 1;2 2.5;1 1.5]; % human speed interval for each section
spd_intv_set = [1 1.5;1 1.5;1 1.5;1 1.5;1 1.5];
dt = 0.5;
samp_t = 0.05;
des_pos = [s_pos,tar_pos]; % human destinations, including his starting position and all target positions
max_iter = 2e4;  % number of points to sample
rand_seed = 30; 

% draw starting position of the human
plot(s_pos(1),s_pos(2),'g^','MarkerSize',10,'MarkerFaceColor','g');

arr_time = zeros(tar_num,1); % time step that robot arrives at a target
traj = [];
for kk = 1:size(des_pos,2)-1
    fld = struct('fld_cor',fld_cor,'obs_set',{obs_set},'g_pos',des_pos(:,kk+1),...
    's_pos',des_pos(:,kk),'obs_margin',obs_margin);
    
%     idx = randi(size(spd_intv_set,1),1);
    spd_intv = spd_intv_set(kk,:);
    step_intv = spd_intv*dt;
    
    display(sprintf('max_iter %d',max_iter))
    display(sprintf('rand_seed %d',rand_seed))
    
    is_demo = 0;
    variant = 'Simple2D'; % this class is used for RRT* for a holonomic vehicle
    
    % find a path using RRT*
    result = pathFinder(fld,step_intv,max_iter,is_demo,rand_seed,variant);
    traj = [traj,result.opt_path];
    arr_time(kk) = size(traj,2);
end

% smooth trajectory
% traj_s = trajSmoothing(traj);
traj_mpc = mpcSmoothing(traj,arr_time,dt,samp_t,spd_intv_set,tar_pos);

% change the sampling time to be samp_t = 0.05 for original trajectory
traj_dense = [];
ratio = dt/samp_t;
for ii = 1:size(traj,2)-1
    pt1 = linspace(traj(1,ii),traj(1,ii+1),ratio+1);
    pt2 = linspace(traj(2,ii),traj(2,ii+1),ratio+1);
    tmp_traj = [pt1(1:end-1);pt2(1:end-1)];
    traj_dense = [traj_dense,tmp_traj];
end
traj_dense = [traj_dense,traj(:,end)]; 

% change the sampling time to be samp_t = 0.05 for smoothed trajectory
% traj_mpc_dense = [];
% ratio = dt/samp_t;
% for ii = 1:size(traj_mpc,2)-1
%     pt1 = linspace(traj_mpc(1,ii),traj_mpc(1,ii+1),ratio+1);
%     pt2 = linspace(traj_mpc(2,ii),traj_mpc(2,ii+1),ratio+1);
%     tmp_traj = [pt1(1:end-1);pt2(1:end-1)];
%     traj_mpc_dense = [traj_mpc_dense,tmp_traj];
% end
% traj_mpc_dense = [traj_mpc_dense,traj(:,end)]; 

% draw trajectory
% hold on
h1 = plot(traj(1,:),traj(2,:),'r');
% h2 = plot(traj_s(1,:),traj_s(2,:),'c');
h3 = plot(traj_mpc(1,:),traj_mpc(2,:),'g');
% h4 = plot(traj_dense(1,:),traj_dense(2,:),'b');
v = sqrt(sum((traj(:,2:end)-traj(:,1:end-1)).^2,1))/dt; % the speed at each time
v3 = sqrt(sum((traj_mpc(:,2:end)-traj_mpc(:,1:end-1)).^2,1))/dt;
%}


%% save the scenario
% scenario = struct('fld_cor',fld_cor,'dt',dt,'s_pos',s_pos,'obs_num',obs_num,...
%     'obs_set',{obs_set},'tar_num',tar_num,'tar_pos',tar_pos,...
%     'obs_margin',obs_margin,'samp_t',samp_t,'traj',traj_dense);
% save('scenario.mat','scenario');
save('workspace.mat');
end

function is_colliding = collisionCheck(shape1, shape2)
% check if two shapes collide

is_colliding = false;

if strcmp(shape1{1},'r')
    w = shape1{2}(1);
    h = shape1{2}(2);
    xmin1 = shape1{3}(1);
    ymin1 = shape1{3}(2);
    xmax1 = shape1{3}(1)+w;
    ymax1 = shape1{3}(2)+h;
elseif strcmp(shape1{1},'c')
    rad1 = shape1{2};
    center1 = shape1{3};
end

if strcmp(shape2{1},'r')
    w = shape2{2}(1);
    h = shape2{2}(2);
    xmin2 = shape2{3}(1);
    ymin2 = shape2{3}(2);
    xmax2 = shape2{3}(1)+w;
    ymax2 = shape2{3}(2)+h;
elseif strcmp(shape2{1},'c')
    rad2 = shape2{2};
    center2 = shape2{3};
end


if strcmp(shape1{1},'r') && strcmp(shape2{1},'r')
    % two rectangles collide iff xmin1<xmax2 && xmax1>xmin2 &&
    % ymin1 < ymax2 && ymax1 > ymin2
    if (xmin1<xmax2 && xmax1>xmin2 && ymin1 < ymax2 && ymax1 > ymin2)
        % colliding happens
        is_colliding = true;
        return
    end
elseif strcmp(shape1{1},'r') && strcmp(shape2{1},'c')
    % the rectangle collide with a circle iff one of two
    % conditions is statisfied:
    % 1. a further point of the circle along x or y
    % direction is inside the rectangle
    % 2. a vertex of rectangle is inside the circle
    
    % test condition 1
    fur_pts = bsxfun(@plus,center2,[-rad2,0,rad2,0;0,rad2,0,-rad2]);
    tmp_cmp = bsxfun(@lt,[xmin1;ymin1],fur_pts) & bsxfun(@gt,[xmax1;ymax1],fur_pts);
    if ~all(sum(tmp_cmp,2)<2)
        % colliding happens
        is_colliding = true;
        return
    end
    
    % test condition 2
    % corner points of the rectangle
    cor_pts = [xmin1,xmin1,xmax1,xmax1;ymin1,ymax1,ymax1,ymin1];
    % distance to circle center
    tmp_cmp = sum((bsxfun(@minus,cor_pts,center2)).^2,1)-rad2^2;
    if ~all(tmp_cmp > 0)
        % colliding happens
        is_colliding = true;
        return
    end
elseif strcmp(shape1{1},'c') && strcmp(shape2{1},'r')
    % the rectangle collide with a circle iff one of two
    % conditions is statisfied:
    % 1. a further point of the circle along x or y
    % direction is inside the rectangle
    % 2. a vertex of rectangle is inside the circle
    
    % test condition 1
    fur_pts = bsxfun(@plus,center1,[-rad1,0,rad1,0;0,rad1,0,-rad1]);
    tmp_cmp = bsxfun(@lt,[xmin2;ymin2],fur_pts) & bsxfun(@gt,[xmax2;ymax2],fur_pts);
    if ~all(sum(tmp_cmp,2)<2)
        % colliding happens
        is_colliding = true;
        return
    end
    
    % test condition 2
    % corner points of the rectangle
    cor_pts = [xmin2,xmin2,xmax2,xmax2;ymin2,ymax2,ymax2,ymin2];
    % distance to circle center
    tmp_cmp = sum((bsxfun(@minus,cor_pts,center1)).^2,1)-rad1^2;
    if ~all(tmp_cmp > 0)
        % colliding happens
        is_colliding = true;
        return
    end
elseif strcmp(shape1{1},'c') && strcmp(shape2{1},'c')
    if norm(center1-center2) <= rad1+rad2
        % colliding happens
        is_colliding = true;
        return
    end
end
end

function traj_s = trajSmoothing(traj)
    max_num = 5;
    len = size(traj,2);
    t = floor(len/max_num);
    traj_s = [];
    
    for ii = 1:t
        if ii < t
            tmp_pts = traj(:,max_num*(ii-1)+1:max_num*ii);            
        else
            tmp_pts = traj(:,max_num*(ii-1)+1:end);
        end
        % polynomial fitting
        p = polyfit(tmp_pts(1,:),tmp_pts(2,:),1);
        fit_y = polyval(p,tmp_pts(1,:));
        traj_s = [traj_s,[tmp_pts(1,:);fit_y]];
    end
end

function traj_mpc = mpcSmoothing(traj,arr_time,dt,samp_t,spd_intv_set,tar_pos)
    traj_len = size(traj,2);
    traj_mpc = traj;
    cnt = 1;
    cur_state = [traj(:,1);0;spd_intv_set(cnt,2)];
%     cur_state = [traj(:,1);max_spd_set(cnt);max_spd_set(cnt)];
    mpc_hor = 10;
    x = sdpvar(4,mpc_hor+1); %[x;y;theta;v]
    u = sdpvar(2,mpc_hor); % [w;a]
    slack = sdpvar(2,mpc_hor);
%     x = sdpvar(4,mpc_hor+1); %[x;y;vx;vy]
%     u = sdpvar(2,mpc_hor); % [ax;ay]
    a_lb = -1;
    a_ub = 1;
    w_lb = -2*pi/3;
    w_ub = 2*pi/3;
%     a_lb = -2;
%     a_ub = 2;
    
    
    for kk = 1:traj_len-mpc_hor
        display(kk)
        v_ub = spd_intv_set(cnt,2);
        v_lb = spd_intv_set(cnt,1);
        
        obj = 0;
        
        if ismember(kk,arr_time)
            cnt = cnt+1;
        end
        
        for ii = 1:mpc_hor
            obj = obj+sum((x(1:2,ii+1)-traj(:,kk+ii)).^2)+0.5*u(1,ii)^2+0.1*u(2,ii)^2;            
        end
        obj = obj++sum(sum(slack.^2));
        constr = [x(:,1) == cur_state, slack>=0];
        
        
        for ii = 1:mpc_hor  
            % double integrator model
%             constr = [constr, x(:,ii+1) == [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1]*x(:,ii) + [0 0;0 0;dt 0;0 dt]*u(:,ii)];
%             constr = [constr, x(3,ii+1)^2+x(4,ii+1)^2 <= max_spd_set(cnt)^2];
%             constr = [constr, a_lb <= u(:,ii) <= a_ub];
            
            % unicycle model
%             constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(4,ii)*[cos(x(3,ii));sin(x(3,ii))]*dt,...
%                 x(3,ii+1) == x(3,ii) + u(1,ii)*dt, x(4,ii+1) == x(4,ii)+u(2,ii)*dt,...
%                 0<=x(4,ii+1)<=max_spd_set(cnt),w_lb<=u(1,ii)<=w_ub,a_lb<=u(2,ii)<=a_ub];
            constr = [constr,x(1:2,ii+1) == x(1:2,ii)+x(4,ii)*[cos(x(3,ii));sin(x(3,ii))]*dt,...
                x(3,ii+1) == x(3,ii) + u(1,ii)*dt, x(4,ii+1) == x(4,ii)+u(2,ii)*dt,...
                spd_intv_set(cnt,1)-slack(1,ii)<=x(4,ii+1)<=spd_intv_set(cnt,2)+slack(2,ii),...
                w_lb<=u(1,ii)<=w_ub,a_lb<=u(2,ii)<=a_ub];
        end
        
        opt = sdpsettings('solver','fmincon','verbose',0);
        sol = optimize(constr,obj,opt);
        cur_state = value(x(:,2));
        if sol.problem ~= 0
            display('mpc goes wrong')
            display(sol.problem)
        end
        if kk < traj_len-mpc_hor
            traj_mpc(:,kk+1) = value(x(1:2,2));
        else
            traj_mpc(:,kk:end) = value(x(1:2,:));
        end
%         display(cur_state)
    end
end