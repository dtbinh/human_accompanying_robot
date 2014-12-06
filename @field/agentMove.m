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
plan_state = inPara.plan_state;
r_state = inPara.r_state;
r_input = inPara.r_input;
k = inPara.k;
hor = inPara.hor;
pre_type = inPara.pre_type;
samp_rate = inPara.samp_rate;
safe_dis = inPara.safe_dis;
mpc_dt = inPara.mpc_dt;
safe_marg = inPara.safe_marg;

%% agents move 
for agentIndex = 1:length(agents)
    agent = agents(agentIndex);
    
    %% human moves
    if strcmp(agent.type,'human')         
        h_next_actions = getNextActionWithFixedHeading(agent.currentPos...
            ,h_tar_wp,agent.currentV,0); % last argument is for zero deg_dev
        
        if k == 1
            tmp_agent_traj = agent.currentPos;
            obv_traj = agent.currentPos(1:2);
        else
            tmp_agent_traj = agent.traj;
        end
        
        % generate observations
        % I thought sometimes the recording of data within 1 second may be
        % less that the sample rate number. So I designed the following
        % code. However, as in endCheck I have already required that the
        % human be considered to be at a target as long as it is within one
        % step range. So in reality the sample number can always equal
        % samp_rate.
        cur_pos = agent.currentPos(1:2);  
        % compute human heading
        if k == 1
            % we assume that human heading at the beginning is known.
            cur_hd = agent.currentPos(3)+h_next_actions(3);
        else
            cur_hd = agent.currentPos(3);
        end
        
        agent = takeNextAction(agent,h_next_actions);
        next_pos = agent.currentPos(1:2);
        t = uint16(norm(cur_pos-next_pos,2)/agent.currentV); % calculate the time for human to move to his next position
        samp_num = double(t*samp_rate*mpc_dt); % get the number of observations of human position
        for ii = 1:samp_num
             % observed human position
            obv_traj = [obv_traj,cur_pos+(next_pos-cur_pos)*ii/samp_num];
        end
        
        % update human position
        agent.traj = [tmp_agent_traj,agent.currentPos];
        agents(agentIndex) = agent;
        
    %% robot predicts and moves
    elseif strcmp(agent.type,'robot')     
        %% estimate human position
        % later should use Donghan's module
        est_state([1,3],k) = obv_traj(:,(k-1)*samp_rate+1); 
        h = agents(1);
        hd = cur_hd;
        est_state([2,4],k) = h.currentV*[cos(hd);sin(hd)];
        
        %%  predict human future path
        % may need to seperately deal with k == 1
        inPara_phj = struct('state',est_state(:,k),'hor',hor,'pre_type',pre_type,...
            'mpc_dt',mpc_dt);
        pre_traj(:,:,k) = predictHumanTraj(agent,inPara_phj);
               
        %% robot path planning
        inPara_pp = struct('pre_traj',pre_traj(:,:,k),'hor',hor,...
            'safe_dis',safe_dis,'mpc_dt',mpc_dt,'h_v',h.currentV,...
            'obs_info',campus.obs_info,'safe_marg',safe_marg);
        outPara_pp = pathPlanner(agent,inPara_pp);
        opt_x = outPara_pp.opt_x;
        opt_u = outPara_pp.opt_u;
        agent.currentPos = opt_x(1:2,2); % robot moves
        agent.currentV = opt_x(3,2); % robot updates its speed
        plan_state(:,:,k) = opt_x;
        
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
% outPara = struct('agents',agents,'obv_traj',obv_traj,'est_state',est_state,...
%     'pre_traj',pre_traj,'plan_state',plan_state,'r_state',fut_x(:,2),'r_input',fut_u(:,1));
outPara = struct('agents',agents,'obv_traj',obv_traj,'est_state',est_state,...
    'pre_traj',pre_traj,'plan_state',plan_state);

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