% 11/24/14
% simulate the movement of each agent
function [outPara] = agentMove(inPara)
%% initialization
% get input arguments
campus = inPara.campus;
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
agentIndex = inPara.agentIndex;
plan_type = inPara.plan_type;
cmft_dis = inPara.cmft_dis;
%% agents move 
% for agentIndex = 1:length(agents)
    agent = agents(agentIndex);
    
    %% human moves
    if strcmp(agent.type,'human')         
        h_next_actions = getNextActionWithFixedHeading(agent.currentPos...
            ,h_tar_wp,agent.currentV,0,mpc_dt); % last argument is for zero deg_dev
        
        if k == 1
            tmp_agent_traj = agent.currentPos;
%             obv_traj = [0;agent.currentPos(1:2)];
        else
            tmp_agent_traj = agent.traj;
        end
        
        % generate observations
        cur_pos = agent.currentPos(1:2);
        %
        % compute human heading
        if k == 1
            % we assume that human heading at the beginning is known.
            cur_hd = agent.currentPos(3)+h_next_actions(3);
        else
            cur_hd = agent.currentPos(3);
        end
        %}
        agent = takeNextAction(agent,h_next_actions);
        next_pos = agent.currentPos(1:2);
        t = norm(cur_pos-next_pos,2)/agent.currentV; % calculate the time for human to move to his next position
        samp_num = double(uint16(t*samp_rate)); % get the number of observations of human position
        for ii = 1:samp_num
             % observed human position
             tmp_t = (k-1)*mpc_dt+(ii-1)/samp_rate;
             obv_traj = [obv_traj,[tmp_t;cur_pos+(next_pos-cur_pos)*(ii-1)/samp_num]];
        end
        
        % update human position
        agent.traj = [tmp_agent_traj,agent.currentPos];
        agents(agentIndex) = agent;
        
    %% robot predicts and moves
    elseif strcmp(agent.type,'robot')  
        guess_u = inPara.guess_u;
        guess_x = inPara.guess_x;
        
        samp_num = inPara.samp_num;
        h = agents(1);
        cur_hd = h.currentPos(3);
        %% estimate human position
        % 
        if strcmp(pre_type,'IMM-UKF')
            [x_est,y_est,x_pos_pre,y_pos_pre] = IMM_UKF_run();
            est_state([1,2],:,k) = x_est((k-1)*samp_num+1:end-1,:)';
            est_state([3,4],:,k) = y_est((k-1)*samp_num+1:end-1,:)';
        elseif strcmp(pre_type,'UKF')
            [x_est,y_est,x_pos_pre,y_pos_pre] = UKF_run();
            est_state([1,2],:,k) = x_est((k-1)*samp_num+1:end-1,:)';
            est_state([3,4],:,k) = y_est((k-1)*samp_num+1:end-1,:)';
        end
        
        %}
        % estimation with no measurement noise
        %{
            est_state([1,3],k) = obv_traj(2:3,(k-1)*samp_rate*mpc_dt+1);
            hd = cur_hd;
            est_state([2,4],k) = h.currentV*[cos(hd);sin(hd)];
        %}
        
        %%  predict human future path
        % prediction by IMM
        %
        if strcmp(pre_type,'IMM-UKF') || strcmp(pre_type,'UKF')
            pre_traj(:,:,k) = [[x_est((k-1)*samp_num+1,1);y_est((k-1)*samp_num+1,1)],[x_pos_pre(k,:);y_pos_pre(k,:)]];
            h_v = [x_est((k-1)*samp_num+1,2);y_est((k-1)*samp_num+1,2)];
%             pre_traj(:,:,k) = [x_pos_pre_imm(:,k)';y_pos_pre_imm(:,k)'];
        % prediction by extrapolation
        elseif strcmp(pre_type,'extpol')
%             inPara_phj = struct('state',est_state(:,k),'hor',hor,'pre_type',pre_type,...
%                 'mpc_dt',mpc_dt);
            inPara_phj = struct('state',[x_est((k-1)*samp_num+1,:)';y_est((k-1)*samp_num+1,:)'],...
                'hor',hor,'pre_type',pre_type,'mpc_dt',mpc_dt);
            pre_traj(:,:,k) = predictHumanTraj(agent,inPara_phj);
            h_v = [x_est((k-1)*samp_num+1,2);y_est((k-1)*samp_num+1,2)];
        elseif strcmp(pre_type,'opt')
            load('sim_res/human_traj.mat','human_traj')
            tmp_len = size(human_traj,2);
            if k+hor <= tmp_len
                pre_traj(:,:,k) = human_traj(1:2,k:k+hor);
            else
                pre_traj(:,:,k) = [human_traj(1:2,k:tmp_len),human_traj(1:2,tmp_len)*ones(1,k+hor-tmp_len)];
            end
            h_v = 1.5;
        end     
%         pos_pre_imm = inPara.pos_pre_imm;
        %}
        %% robot path planning
        %
        if strcmp(plan_type,'MPC')         
            %%% a possible bug: h_v here is a vector while in
            %%% pathPlanner.m, it is treated as a scalar.
            inPara_pp = struct('pre_traj',pre_traj(:,:,k),'hor',hor,...
                'safe_dis',safe_dis,'mpc_dt',mpc_dt,'h_v',h_v,...
                'obs_info',campus.obs_info,'safe_marg',safe_marg,...
                'guess_u',guess_u,'guess_x',guess_x,'cmft_dis',cmft_dis);

            % assume accurate estimation and prediction
            %{
            load('human_traj.mat','human_traj');
            inPara_pp = struct('pre_traj',human_traj(1:2,k:k+hor),'hor',hor,...
                'safe_dis',safe_dis,'mpc_dt',mpc_dt,'h_v',1.5,...
                'obs_info',{campus.obs_info},'safe_marg',safe_marg,...
                'guess_u',guess_u,'guess_x',guess_x);
            %}
            outPara_pp = pathPlanner(agent,inPara_pp);
            opt_x = outPara_pp.opt_x;
            opt_u = outPara_pp.opt_u;
            agent.currentPos = opt_x(1:3,2); % update robot position
            agent.currentV = opt_x(4,2); % update robot speed
            r_state(:,k+1) = opt_x(:,2);
            r_input(:,k) = opt_u(:,1);
            plan_state(:,:,k) = opt_x;
            
            guess_u = [opt_u(:,2:end),zeros(size(opt_u,1),1)];
            guess_x = [opt_x(:,2:end),opt_x(:,end)];
        elseif strcmp(plan_type,'greedy1') || strcmp(plan_type,'greedy0')
%             inPara_pp = struct('pre_traj',pre_traj(:,:,k),'hor',hor,...
%                 'safe_dis',safe_dis,'mpc_dt',mpc_dt,'h_v',...
%                 [x_est((k-1)*samp_num+1,2);y_est((k-1)*samp_num+1,2)],'obs_info',campus.obs_info,...
%                 'safe_marg',safe_marg,'plan_type',plan_type);

            % assume accurate estimation and prediction
%             load('human_traj.mat','human_traj');
            inPara_pp = struct('pre_traj',pre_traj(:,1,k),'hor',2,...
                'safe_dis',safe_dis,'mpc_dt',mpc_dt,'h_v',h_v,...
                'obs_info',campus.obs_info,'safe_marg',safe_marg,...
                'guess_u',guess_u,'guess_x',guess_x,'cmft_dis',cmft_dis);
            outPara_pp = pathPlannerGreedy(agent,inPara_pp);
            opt_x = outPara_pp.opt_x;
            opt_u = outPara_pp.opt_u;
            agent.currentPos = opt_x(1:3,2); % robot moves
            agent.currentV = opt_x(4,2); % robot updates its speed
            r_state(:,k+1) = opt_x(:,2);
            r_input(:,k) = opt_u(:,1);
            plan_state(:,1:size(opt_x,2),k) = opt_x;
        end
        %}
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
% end
%% define output arguments
% outPara = struct('agents',agents,'obv_traj',obv_traj,'est_state',est_state,...
%     'pre_traj',pre_traj,'plan_state',plan_state,'r_state',fut_x(:,2),'r_input',fut_u(:,1));
outPara = struct('agents',agents,'obv_traj',obv_traj,'est_state',est_state,...
    'pre_traj',pre_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
    'samp_num',samp_num);
if exist('guess_u', 'var')
    outPara.guess_u = guess_u;
end  
if exist('guess_x', 'var')
    outPara.guess_x = guess_x;
end 
end

function next_act = getNextActionWithFixedHeading(a_pos,t_pos,v,deg_dev,mpc_dt)
% calculate the actions for agent to move from a_pos to target position
% t_pos with velocity v
vec = t_pos - a_pos(1:2);
heading = calAngle(vec)+deg_dev;
% next_acts = zeros(3); %[dx,dy,d_heading]. only the first action changes the agent's heading
tmp_a_pos = a_pos;
dx = v*cos(heading)*mpc_dt;
dy = v*sin(heading)*mpc_dt;
d_hd = heading - tmp_a_pos(3);
next_act = [dx;dy;d_hd];
end    