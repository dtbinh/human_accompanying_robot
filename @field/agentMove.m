% 11/24/14
% simulate the movement of each agent
function [outPara] = agentMove(campus,inPara)
%% initialization
% get input arguments
% campus = inPara.campus;
agents = inPara.agents;
% way_pts = inPara.way_pts;
h_tar_wp = inPara.h_tar_wp;
obv_traj = inPara.obv_traj;
est_state = inPara.est_state;
pre_traj = inPara.pre_traj;
plan_path = inPara.plan_path;
pre_type = inPara.pre_type;
r_path = inPara.r_path;
k = inPara.k;
hor = inPara.hor;

%% agents move 
for agentIndex = 1%1:length(agents)
    agent = agents(agentIndex);
    
    %% human moves
    if strcmp(agent.type,'human')         
        h_next_actions = getNextActionWithFixedHeading(agent.currentPos...
            ,h_tar_wp,agent.maxV,0); % last argument is for zero deg_dev
        % Take action
        if k == 1
            tmp_agent_traj = agent.currentPos;
        else
            tmp_agent_traj = agent.traj;
        end
                
        agent = takeNextAction(agent,h_next_actions);

        % Save results
        obv_traj(:,k) = agent.currentPos(1:2); % observed human position, later should add noise
        agent.traj = [tmp_agent_traj,agent.currentPos];
        agents(agentIndex) = agent;
        
    %% robot predicts and moves
    elseif strcmp(agent.type,'robot')     
        %% estimate human position
        % later should use Donghan's module
        est_state([1,3],k) = obv_traj(:,k); 
        
        %%  predict human future path
        % may need to seperately deal with k == 1
        inPara_phj = struct('state',est_state(:,k),'hor',hor,'pre_type',pre_type);
        pre_traj(:,:,k) = predictHumanTraj(agent,inPara_phj);
               
        %% robot path planning
        inPara_pp = struct('pre_traj',pre_traj(:,:,k));
        plan_path(2,:,k) = pathPlanner(agent,inPara_pp);
        r_path = plan_path(2,1,k);
        agent = takeNextAction(agent, r_path);

        if k == 1
            tmp_agent_traj = agent.currentPos;
        else
            tmp_agent_traj = agent.traj;
        end

        agent.traj = [tmp_agent_traj,agent.currentPos];
        agents(agentIndex) = agent;

    else
        error('Invalid agent type for planning path')
    end
end
%% define output arguments
outPara = struct('agents',agents,'obv_traj',obv_traj,'est_state',est_state,...
    'pre_traj',pre_traj,'plan_path',plan_path,'r_path',r_path);
end

function next_act = getNextActionWithFixedHeading(a_pos,t_pos,v,deg_dev)
% calculate the actions for agent to move from a_pos to target position
% t_pos with velocity v
vec = t_pos - a_pos(1:2);
heading = calAngle(vec)+deg_dev;
% next_acts = zeros(3); %[dx,dy,d_heading]. only the first action changes the agent's heading
tmp_a_pos = a_pos;
dx = v*cos(heading);
dy = v*sin(heading);
d_hd = heading - tmp_a_pos(3);
next_act = [dx;dy;d_hd];
end
    