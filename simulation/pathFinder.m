function result = pathFinder(fld,step_intv,max_iter,is_demo,rand_seed,variant)
% use RRT* to find path
%%% Configuration block
conf = struct;
conf.delta_goal_point = 1;          % Radius of goal point
conf.delta_near = 3;%max_step;              % Radius for neighboring nodes % 1.5
% conf.max_step = max_step;                % Maximum position change when we add a new node to the tree
conf.step_intv = step_intv;           % speed interval, representing the min and max speed
conf.bin_size = 7;
conf.demo = is_demo;

% frame_set(max_iter+1) = struct('cdata',[],'colormap',[]); % save plot frame for animation

hold on
% plot(fld.s_pos(1),fld.s_pos(2),'g^','MarkerSize',10,'MarkerFaceColor','g');

MAX_NODES   = max_iter;
MAX_ITER    = max_iter;
RAND_SEED   = rand_seed;
MAP = fld;
CONF = conf;

addpath(genpath(pwd));

% ALGORITHM = 'RRTS';

% call the class simple2D here
problem = eval([variant '(RAND_SEED, MAX_NODES, MAP, CONF);']);
%problem = RedundantManipulator(RAND_SEED, MAX_NODES, MAP, CONF);

% disp(ALGORITHM);

% draw starting point
% if is_demo == 1
%     plot(fld.s_pos(1),fld.s_pos(2),'ro')
%     hold on
% %     frame_set(1) = getframe(gcf);
% end

% % starting timer
tic;
for ind = 1:MAX_ITER-1
    new_node = problem.sample();
    nearest_node_ind = problem.nearest(new_node);
    new_node = problem.steer(nearest_node_ind, new_node);   % if new node is very distant from the nearest node we go from the nearest node in the direction of a new node
    
    % plot newly generated new node
    if is_demo == 1
        plot(new_node(1),new_node(2),'ro')
        hold on
    end
    
    if(~problem.collision_checker(new_node, nearest_node_ind))
        neighbors = problem.neighbors(new_node, nearest_node_ind);
        min_node_ind = problem.chooseParent(neighbors, nearest_node_ind, new_node);
        new_node_ind = problem.insert_node(min_node_ind, new_node);
        
        % draw the edge between min_node and new_node
        if is_demo == 1
            hdl = plot ([problem.tree(1,min_node_ind),problem.tree(1,new_node_ind)],[problem.tree(2,min_node_ind),problem.tree(2,new_node_ind)]);
            edge_pair = sprintf('%da%d',min_node_ind,new_node_ind);
            problem.handle(edge_pair) = hdl;
            drawnow
        end        
        problem.rewire(new_node_ind, neighbors, min_node_ind);
        
    end
    
    if is_demo == 1
        if(mod(ind, 1000) == 0)
            disp([num2str(ind) ' iterations ' num2str(problem.nodes_added-1) ' nodes in ' num2str(toc) ' rewired ' num2str(problem.num_rewired)]);
        end        
    end
end

[opt_path,is_reach] = problem.plot();


result = struct;
result.is_reach = is_reach;
result.opt_path = opt_path;
end