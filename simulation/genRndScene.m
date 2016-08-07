function genRndScene()
%%% randomly generate scenarios for simulation
clc;
clear all;

% basic setup
fld_cor = [0,0,100,100]; % size of the field [xmin,ymin,xmax,ymax]
s_pos = [5;5]; % starting position of the human

%% define obstacle
% define basic obstacle shapes
% {'r',[w,h]}: rectangle with horizontal size a (width) and vertical size b
% (height).
% {'c',a}: circle with radius a

shape_set = {{'r',[5,20]},{'r',[20,5]},{'r',[5,30]},{'c',5}}; 

% generate obstacle positions
obs_num = 5; % specifiy total number of obstacles
obs_set = cell(obs_num,1);
count = 0; % counter for obstacle

while (count < obs_num)
    shape_idx = randi(length(shape_set),1); % randomly choose a shape    
    tmp_shape = shape_set{shape_idx};
    
    tmp_pos = [randi([fld_cor(1),fld_cor(3)],1);randi([fld_cor(2),fld_cor(4)],1)]; % find its position
    if strcmp(tmp_shape{1},'c')        
    % for circle, the generated position is its center
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

obs_margin = 3; % the target should be away from the boundaries of obstacles by this margin

% specify target number
tar_num = 5; 
tar_cnt = 1;
tar_pos = zeros(2,tar_num);

while (tar_cnt <= tar_num)
    re_gen = false; % if true, the generated point is inside an obstacle and new regeneration
    
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
            w = tmp_shape{2}(1);
            h = tmp_shape{2}(2);
            ll_cor = tmp_shape{3}; % lower-left coordinate
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
% RRT* 
hold on
addpath('rrtstar')
% plot(s_pos(1),s_pos(2),'g^','MarkerSize',10,'MarkerFaceColor','g');

max_spd = 1; % human speed
dt = 0.5;
max_step = max_spd*dt;
des_pos = [s_pos,tar_pos]; % human destinations, including his starting position and all target positions
max_iter = 9e3;  % number of points to sample
rand_seed = 30; 

fld = struct('fld_cor',fld_cor,'obs_set',{obs_set},'g_pos',zeros(2,1),...
    's_pos',zeros(2,1),'obs_margin',obs_margin);
traj = [];
for kk = 1:size(des_pos,2)-1
    fld.s_pos = des_pos(:,kk);
    fld.g_pos = des_pos(:,kk+1);
    
    display(sprintf('max_iter %d',max_iter))
    display(sprintf('rand_seed %d',rand_seed))
    
    is_demo = 0;
    variant = 'Simple2D'; % this class is used for RRT* for a holonomic vehicle
    
    result = pathFinder(fld,max_step,max_iter,is_demo,rand_seed,variant);
    traj = [traj,result.opt_path];
end
% draw trajectory
plot(traj(1,:),traj(2,:),'r')

%}

% MPC
%{
% use MPC to find feasible path
h_max_spd = 1; % human speed
dt = 0.5;
des_pos = [s_pos,tar_pos]; % human destinations, including his starting position and all target positions
hor = 50;
rect_num = 0;

% find the number of rectangles
for jj = 1:obs_num 
    tmp_shape = obs_set{jj};
    if strcmp(tmp_shape{1},'r')
        rect_num = rect_num+1;
    end
end

x = sdpvar(2,hor+1); % [x;y] coordinate of human
u = sdpvar(2,hor); % [vx;vy]
lambda = sdpvar(rect_num,4,hor); % dummy variables

% constraints
for ii = 1:hor
    % kinematics
    constr0 = [x(:,ii+1) == x(:,ii)+u(:,ii)*dt];
    
    % collision avoidance
    rect_cnt = 1; % counter for rectangular obstacles
    
    % x > xmax+obs_margin || x < xmin-obs_margin || y > ymax+obs_margin || y < ymin-obs_margin  
    for jj = 1:obs_num %
        tmp_shape = obs_set{jj};
        if strcmp(tmp_shape{1},'c')
            % for circle, the generated position is its center
            radius = tmp_shape{2};
            constr0 = [constr0, sum((x(:,ii+1)-tmp_shape{3}).^2) <= (radius+obs_margin)^2];
        elseif strcmp(tmp_shape{1},'r')
            % for rectangle, the generated position is its lower-left coordinate
            w = tmp_shape{2}(1);
            h = tmp_shape{2}(2);
            ll_cor = tmp_shape{3}; % lower-left coordinate
            xmin = ll_cor(1);
            xmax = ll_cor(1)+w;
            ymin = ll_cor(2);
            ymax = ll_cor(2)+h;
            constr0 = [constr0, lambda(rect_cnt,:,ii)*[xmax+obs_margin-x(1,ii+1);...
                x(1,ii+1)-xmin+obs_margin;ymax+obs_margin-x(2,ii+1);...
                x(2,ii+1)-ymin+obs_margin] <= 0];
            rect_cnt = rect_cnt+1;
        end
    end
    
    constr0 = [constr0,sum(squeeze(lambda(:,:,ii)),2) == ones(size(lambda,1),1)];
end

% bounds on state and input
constr0 = [constr0, repmat(fld_cor(1:2)',1,size(x,2)) <= x <= repmat(fld_cor(3:4)',1,size(x,2))];
constr0 = [constr0, sum(u.^2,1) <= h_max_spd^2];
constr0 = [constr0, lambda >= 0];

% solve for each pair of start and goal position
for ii = 1%:size(des_pos,2)-1  
    % initial condition
    constr = [constr0,x(:,1) == des_pos(:,ii)];
    
    % initial guess
    dx = (des_pos(:,ii+1)-des_pos(:,ii))/hor;
    guess_x = [des_pos(1,ii):dx(1):des_pos(1,ii+1);des_pos(2,ii):dx(2):des_pos(2,ii+1)];
    guess_u = guess_x(:,2:end)-guess_x(:,1:end-1);
    
    % obj
    obj = sum(sum(([x(:,2:end),des_pos(:,ii+1)]-x(:,1:end)).^2,1))+...
        sum(sum(u.^2,1));

    % solve MPC
    % initial guess

    if ~isempty(guess_x)
        assign(x,guess_x); % assign initial solution
        assign(u,guess_u); % assign initial input
    end

%     opt = sdpsettings('solver','fmincon','fmincon.Algorithm','interior-point','usex0',1,'debug',0,'verbose',1);    
    opt = sdpsettings('solver','ipopt','usex0',1,'debug',1,'verbose',1);    
    
    % solve nlp
    tic
    sol = optimize(constr,obj,opt);
    elapse_time = toc
    
    optx = value(x);
end

plot(optx(1,:),optx(2,:))

%}

% find path to bypass rectangle
%{
h_spd = 1; % human speed
dt = 0.5;
des_pos = [s_pos;tar_pos]; % human destinations, including his starting position and all target positions
waypts = s_pos;% waypoints
for ii = 1:size(des_pos,2)-1    
    tmp_s = des_pos(:,ii); % start
    tmp_g = des_pos(:,ii+1); % end
    cur_waypt = tmp_s;
    ang = calAngle(tmp_g-tmp_s);
    
    % find the next way point until human arrives at the goal
    while (norm(cur_waypt-tmp_g) >= v*dt)
        next_waypt = cur_waypt+v*[cos(ang);sin(ang)]*dt;
        for jj = 1:obs_num % check if the generated point is in an obstacle
            tmp_shape = obs_set{jj};
            if strcmp(tmp_shape{1},'c')
                % for circle, the generated position is its center
                radius = tmp_shape{2};
                if norm(next_waypt-tmp_shape{3}) <= radius+obs_margin
                    tmp_waypts = bypassCircle;
                    
                end
            elseif strcmp(tmp_shape{1},'r')
                % for rectangle, the generated position is its lower-left coordinate
                w = tmp_shape{2}(1);
                h = tmp_shape{2}(2);
                ll_cor = tmp_shape{3}; % lower-left coordinate
                if (next_waypt(1) >= ll_cor(1)-obs_margin) && (next_waypt(2) >= ll_cor(2)-obs_margin) && ...
                        (next_waypt(1) <= ll_cor(3)+obs_margin) && (next_waypt(2) <= ll_cor(4)+obs_margin)
                    re_gen = true;
                    break
                end
            end
        end
    end
end
%}
end

% function [way_pt] = findWayPoint()
% % if a line connecting two waypoints collide with an obstacle, define a new
% % way point to by pass the obstacle
% 
% 
% end
% 
% function [is_collide,col_obs,col_pt] = collisionChecker()
% % check whether the staightline between two waypoints intersect with any
% % obstacles. If yes, find the colliding obstacle and the point of contact
% 
% 
% end

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