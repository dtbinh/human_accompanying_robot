% 11/24/2014
% this file is for the ME 290J course project
clc 
clear % clear global variables
close all

%% Setup
scale = 1/3; % scale the size of the field
set(0,'DefaultFigureWindowStyle','docked');
%%% define agents %%%
% Human agent 1
h = agent('human');
h.currentPos = [30;10;0]*scale;%[290;30;0]; % [x y heading]
% h.maxV = 1.5;
h.currentV = 1.5;

% Robot agent 1
r = agent('robot');
r.currentPos = [20;10;0]*scale;%[310;30;0]; %[23.5;0.5;0];
r.currentV = 1;
r.a_lb = -3; 
r.a_ub = 1;
r.w_lb = -pi/4;
r.w_ub = pi/4;

%%% Set field %%%
xLength = 300*scale; 
yLength = 300*scale; 
xMin = 0;
yMin = 0;

%%% set way points
tar_pos = [90,230,210,250,20;60,30,130,230,260]*scale; % target positions
step_size = 1;
% manually pick the way pts for simulated human
inPara_gwp = struct('tar_pos',tar_pos,'scale',scale,'type','h');
h_way_pts = getWayPts(inPara_gwp);
% apply different human speed between different way points.
% h_v = [2,3,1,1,1,1,1,1,1,1,3,1.5,2,3,2,1.5,4];
% apply different acceleration for human speed
% h_acl = -h.maxA+2*h.maxA*rand(1,300);
%%%

% define vertices of obstacles. Here use round obstacles
c_set = [100,65,0,220*scale;33,100,65,40*scale];
r_set = [20,20,20,10*scale];
theta_set = {{pi/2:pi/8:3*pi/2;pi:pi/8:2*pi;-pi/2:pi/8:pi/2;0:pi/8:2*pi}};
inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta_set',theta_set,'type','obs');
obs_pos = getWayPts(inPara_gwp);

campus = field([xMin xMin+xLength/step_size yMin yMin+yLength/step_size],tar_pos);
campus.agentNum = 2;
campus.obs_info = [c_set;r_set]; % gives the center and radius of each obstacle
%% Simulation
% simulation parameters
kf = 500; % simulation length (/s)
agents = [h r];
hor = 5; % MPC horizon 
pre_type = 'IMM';%'extpol'; % 'extpol','IMM'. specify the method for predicting human motion
plan_type = 'MPC'; % 'MPC','greedy1','greedy0'. specify the method for robot controller
samp_rate = 20; % sampling rate (/Hz)
safe_dis = 2; %safe distance between human and robot
safe_marg = 2; % safety margin between human the the obstacle
mpc_dt = 0.5; % sampling time for model discretization used in MPC

% initialize variables
obv_traj = zeros(3,0); % observed human trajectory; first row denotes the time [t,x,y]
est_state = zeros(4,mpc_dt*samp_rate,kf); % estimated human states for every second [x,vx,y,vy];
pre_traj = zeros(2,hor+1,kf); % current and predicted future human trajectory [x,y]
plan_state = zeros(4,hor+1,kf); % robot's current and planned future state [x,y,v]
r_state = zeros(4,kf); % robot's actual state [x,y,theta,v]
r_input = zeros(2,kf); % robot's actual control input [w,a]
wp_cnt = 1; % the waypoint that the human is heading for
h_tar_wp = h_way_pts(:,wp_cnt); % the way point that the human is heading for

for k = 1:kf
    display(k)
    
    %% change human speed
    %{
    tmp_v = agents(1).currentV+h_acl(k)*mpc_dt;
    if tmp_v <= 0
        tmp_v = 0;
    else
        agents(1).currentV = tmp_v;
    end
    %}
    
    %% waypoint check
    % check if the human needs to choose the next way point
    inPara_ec = struct('h',agents(1),'h_way_pts',h_way_pts,'wp_cnt',wp_cnt,'mpc_dt',mpc_dt);
    [outPara_ec] = endCheck(inPara_ec);

    game_end = outPara_ec.game_end;
    arr_wp_flag = outPara_ec.arr_wp_flag; % the human has arrived at a way point
    
    if game_end == 1
        sprintf('total search time is %d',k-1)
        break
    end
    
    if arr_wp_flag == 1
        wp_cnt = wp_cnt+1;
        h_tar_wp = h_way_pts(:,wp_cnt);
%         agents(1).currentV = h_v(wp_cnt);
    end
    
    %% both agents move
    % human moves to a way point and robot predicts and plans its own path
    % human moves
    agentIndex = 1;
    inPara_ams = struct('campus',campus,'agents',agents,'h_tar_wp',h_tar_wp,...
        'obv_traj',obv_traj,'est_state',est_state,...
        'pre_traj',pre_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
        'k',k,'hor',hor,'pre_type',pre_type,'samp_rate',samp_rate,...
        'safe_dis',safe_dis,'mpc_dt',mpc_dt,'safe_marg',safe_marg,...
        'agentIndex',agentIndex,'plan_type',plan_type);
    [outPara_ams] = agentMove(inPara_ams);
    agents = outPara_ams.agents;
    obv_traj = outPara_ams.obv_traj;
    samp_num = outPara_ams.samp_num;
    
    % robot moves
    %
    agentIndex = 2;
%     load('obv_traj3_w_time.mat')% Load Tracjectory of Human
    obv_traj1=obv_traj';
    parameter;  % parameter for IMM
    inPara_ams = struct('campus',campus,'agents',agents,'h_tar_wp',h_tar_wp,...
        'obv_traj',obv_traj1','est_state',est_state,...
        'pre_traj',pre_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
        'k',k,'hor',hor,'pre_type',pre_type,'samp_rate',samp_rate,...
        'safe_dis',safe_dis,'mpc_dt',mpc_dt,'safe_marg',safe_marg,...
        'agentIndex',agentIndex,'plan_type',plan_type,'samp_num',samp_num);
    [outPara_ams] = agentMove(inPara_ams);
    agents = outPara_ams.agents;
    est_state = outPara_ams.est_state;
    pre_traj = outPara_ams.pre_traj;
    plan_state = outPara_ams.plan_state; 
    r_state = outPara_ams.r_state;
    r_input = outPara_ams.r_input;
    %}
    %% plot trajectories
    % Plot future agent positions
    
    % plot specifications
    color_agent = {'r','g','r','g'};
    marker_agent = {'o','^','*','d'};
    line_agent = {'-','-','-','-'};
    orange = [1 204/255 0];
    color_target = {'m','b',orange};
    figure;
    hold on

    % draw targets
    for jj = 1:campus.targetNum
        h = plot(campus.targetPos(1,jj),campus.targetPos(2,jj),'MarkerSize',11);
        set(h,'Marker','p');
    end
    
    % draw obstacles
    for jj = 1:size(obs_pos)
        tmp_pos = obs_pos{jj};
        fill(tmp_pos(1,:),tmp_pos(2,:),'b');
    end
    
    xlim([0,campus.endpoints(2)]);
    ylim([0,campus.endpoints(4)]);
    
    % draw agent trajectory
    for ii = 1:length(agents)
        tmp_agent = agents(ii);
        h1 = plot(tmp_agent.traj(1,:),tmp_agent.traj(2,:),'markers',1);
        set(h1,'MarkerFaceColor',color_agent{ii});
        set(h1,'MarkerEdgeColor',color_agent{ii});
        set(h1,'Color',color_agent{ii});
        set(h1,'LineStyle',line_agent{ii});
        set(h1,'Marker',marker_agent{ii});
        h2 = plot(tmp_agent.currentPos(1),tmp_agent.currentPos(2),color_agent{ii},'markers',3);
        set(h2,'MarkerFaceColor',color_agent{ii});
        set(h2,'MarkerEdgeColor',color_agent{ii});
        set(h2,'Color',color_agent{ii});
        set(h2,'LineStyle',line_agent{ii});
        set(h2,'Marker',marker_agent{ii});
    end
    % predicted human positions
    h3 = plot(pre_traj(1,:,k),pre_traj(2,:,k),color_agent{3},'markers',1);
    set(h3,'MarkerFaceColor',color_agent{3});
    set(h3,'MarkerEdgeColor',color_agent{3});
    set(h3,'Color',color_agent{3});
    set(h3,'LineStyle',line_agent{3});
    set(h3,'Marker',marker_agent{3});
    
    % planned robot trajectory
    h4 = plot(plan_state(1,:,k),plan_state(2,:,k),color_agent{4},'markers',1);
    set(h4,'MarkerFaceColor',color_agent{4});
    set(h4,'MarkerEdgeColor',color_agent{4});
    set(h4,'Color',color_agent{4});
    set(h4,'LineStyle',line_agent{4});
    set(h4,'Marker',marker_agent{4});
    grid minor
    set(gca,'MinorGridLineStyle','-','XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5])
    axis equal
    xlim([0,xLength]);ylim([0,yLength]);
    
    % close plots when there are too many plots
    h5 = gcf;
    if h5 > 50 
        close all;
    end
    %}
end

%% save simulation result
%{
% save data
% if the data is a decimal, replace the '.' with 'p'
str_safe_dis = strrep(num2str(safe_dis),'.','p');
str_safe_marg = strrep(num2str(safe_marg),'.','p');
str_h_v = strrep(num2str(agents(1).currentV),'.','p');

str_t = datestr(clock,'dd-mmm-yyyy_HHMMSS');
folder_path = ('.\sim_res');
data_name = sprintf('sim_traj_%s_%s_%s_%s_%s_%s.mat',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name = fullfile (folder_path,data_name);
save(file_name,'obv_traj','est_state','pre_traj','plan_state','r_state','r_input');

% save plot
folder_path = ('.\sim_res');
fig_name = sprintf('sim_traj_%s_%s_%s_%s_%s_%s.fig',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name = fullfile (folder_path,fig_name);
h = gcf;
saveas (h,file_name);

% convert .fig to .pdf
fig_name2 = sprintf('sim_traj_%s_%s_%s_%s_%s_%s.pdf',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name2 = fullfile (folder_path,fig_name2);
fig2Pdf(file_name2,300,h)
%}