function genRndScene()
% randomly generate scenarios

% basic setup
fld_cor = [0,0,100,100]; % [xmin,ymin,xmax,ymax]
s_pos = [5;5]; % starting position

%% define obstacle
% define basic obstacle shapes
% {'r',[w,h]}: rectangle with horizontal size a (width) and vertical size b
% (height).
% {'c',a}: circle with radius a

shape_set = {{'r',[5,20]},{'r',[20,5]},{'r',[5,30]},{'c',5}}; 

% generate obstacle positions
obs_num = 6; % specifiy total number of obstacles
obs_set = cell(obs_num,1);
count = 1; % index of obstacle
while (count <= obs_num)
    shape_idx = randi(length(shape_set),1); % randomly choose a shape    
    tmp_shape = shape_set{shape_idx};
    
    tmp_pos = [randi([fld_cor(1),fld_cor(3)],1);randi([fld_cor(2),fld_cor(4)],1)]; % find its position
    if strcmp(tmp_shape{1},'c')        
    % for circle, the generated position is its center
        radius = tmp_shape{2};
        
        % keep this obstacle if it falls within the field. otherwise,
        % generate a new one
        if (tmp_pos(1)+radius <= fld_cor(3)) && (tmp_pos(1)-radius >= fld_cor(1)) ...
                && (tmp_pos(2)+radius <= fld_cor(4)) && (tmp_pos(2)-radius >= fld_cor(2))
            obs_set{count} = {shape_set{shape_idx}{1},shape_set{shape_idx}{2},tmp_pos};
            count = count + 1;
        end
    elseif strcmp(tmp_shape{1},'r')
        % for rectangle, the generated position is its lower-left coordinate
        w = tmp_shape{2}(1);
        h = tmp_shape{2}(2);
        
        % keep this obstacle if it falls within the field. otherwise,
        % generate a new one
        if (tmp_pos(1)+w <= fld_cor(3)) && (tmp_pos(2)+h <= fld_cor(4))
            obs_set{count} = {shape_set{shape_idx}{1},shape_set{shape_idx}{2},tmp_pos};
            count = count + 1;
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
        plot(center(1)+radius*cos(delta_ang),center(2)+radius*sin(delta_ang));
        
    elseif strcmp(tmp_shape{1},'r')
        % for rectangle, the generated position is its lower-left coordinate
        w = tmp_shape{2}(1);
        h = tmp_shape{2}(2);
        ll_cor = tmp_shape{3}; % lower-left coordinate
        all_cor = [ll_cor,ll_cor+[w;0],ll_cor+[w;h],ll_cor+[0;h]];
        plot([all_cor(1,:),all_cor(1,1)],[all_cor(2,:),all_cor(2,1)]);
    end
end

% draw targets
for ii = 1:size(tar_pos,2)
    plot(tar_pos(1,ii),tar_pos(2,ii),'r*')
end

%% find feasible path between targets
%%% not finished yet. resume from here. Probably use recursion to generate
%%% path, or use A*

% find path to bypass rectangle
%
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

function [way_pt] = findWayPoint()
% if a line connecting two waypoints collide with an obstacle, define a new
% way point to by pass the obstacle


end

function [is_collide,col_obs,col_pt] = collisionChecker()
% check whether the staightline between two waypoints intersect with any
% obstacles. If yes, find the colliding obstacle and the point of contact


end