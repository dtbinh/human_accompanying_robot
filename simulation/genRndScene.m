function genRndScene()
%%% randomly generate scenarios for simulation
% Chang Liu 8/6/16
% parameters to tune:
% s_pos: starting position of the human
% shape_set: the set of obstacles shapes
% obs_num: total number of obstacles
% tar_num: total number of targets
% max_spd_set: the maximum speed of the human between a pair of obstacles

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

obs_margin = 3; % the target should be away from the boundaries of obstacles by this margin

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
% this code usese RRT* to find a path
hold on
addpath('rrtstar')

max_spd_set = [1,2,1,3,2]; % human speed for each section
dt = 0.5;
samp_t = 0.05;
des_pos = [s_pos,tar_pos]; % human destinations, including his starting position and all target positions
max_iter = 9e3;  % number of points to sample
rand_seed = 30; 

% draw starting position of the human
plot(s_pos(1),s_pos(2),'g^','MarkerSize',10,'MarkerFaceColor','g');

traj = [];
for kk = 1:size(des_pos,2)-1
    fld = struct('fld_cor',fld_cor,'obs_set',{obs_set},'g_pos',des_pos(:,kk+1),...
    's_pos',des_pos(:,kk),'obs_margin',obs_margin);
    
    max_spd = max_spd_set(kk);
    max_step = max_spd*dt;
    
    display(sprintf('max_iter %d',max_iter))
    display(sprintf('rand_seed %d',rand_seed))
    
    is_demo = 0;
    variant = 'Simple2D'; % this class is used for RRT* for a holonomic vehicle
    
    % find a path using RRT*
    result = pathFinder(fld,max_step,max_iter,is_demo,rand_seed,variant);
    traj = [traj,result.opt_path];   
end

% smooth trajectory
traj_s = trajSmoothing(traj);
assignin('base','traj',traj);
assignin('base','traj_s',traj_s);

% change the sampling time to be samp_t = 0.05
traj_dense = [];
ratio = dt/samp_t;
for ii = 1:size(traj_s,2)-1
    pt1 = linspace(traj_s(1,ii),traj_s(1,ii+1),ratio);
    pt2 = linspace(traj_s(2,ii),traj_s(2,ii+1),ratio);
    tmp_traj = [pt1(1:end-1);pt2(1:end-1)];
    traj_dense = [traj_dense,tmp_traj];
end

% draw trajectory
hold on
plot(traj(1,:),traj(2,:),'r')
plot(traj_dense(1,:),traj_dense(2,:),'b')
%}

%% save the scenario
% scenario = struct('fld_cor',fld_cor,'dt',dt,'s_pos',s_pos,'obs_num',obs_num,...
%     'obs_set',{obs_set},'tar_num',tar_num,'tar_pos',tar_pos,...
%     'obs_margin',obs_margin,'samp_t',samp_t,'traj',traj_dense);
% save('scenario.mat','scenario');
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
    max_num = 20;
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
        p = polyfit(tmp_pts(1,:),tmp_pts(2,:),3);
        fit_y = polyval(p,tmp_pts(1,:));
        traj_s = [traj_s,[tmp_pts(1,:);fit_y]];
    end
end