% 11/24/2014
% this file is for the ME 290J course project
clc 
clear % clear global variables
close all

%% Setup
addpath('/Users/changliu/Documents/MATLAB/studentSnopt')
% include IPOPT in YALMIP
% addpath('D:\Program Files\MATLAB\2013a_crack\IPOPT3.11.8');
% addpath('D:\Chang Liu\ipopt');
addpath('/Users/changliu/Documents/MATLAB/Ipopt-3.11.8-linux64mac64win32win64-matlabmexfiles')

scale = 1/3; % scale the size of the field
set(0,'DefaultFigureWindowStyle','docked');
%%% define agents %%%
% Human agent 1
h = agent('human');
h.currentPos = [5;10/3;0]; %[10;10/3;0] [x y heading]
% h.maxV = 1.5;
h.currentV = 1.5;

% Robot agent 1
r = agent('robot');
r.currentPos = [5/3;10/3;0];%[20/3;10/3;0]; 
r.currentV = 1;
r.a_lb = -3; 
r.a_ub = 1;
r.w_lb = -pi/6;
r.w_ub = pi/6;

%%% Set field %%%
x_dis = 16; % displacement on x direction. x_dis>0 means the original of x moves right
y_dis = 18; % displacement on y direction. x_dis>0 means the top of y moves downward (e.g., y_dis = 10 means, if the yLength = 100, the maximum y coordiante is 90)
xLength = 100 - x_dis;
yLength = 100 - y_dis; 
xMin = 0;
yMin = 0;

%%% set way points
tar_pos = [30-x_dis,230/3-x_dis,70-x_dis,250/3-x_dis,44-x_dis,44-x_dis;20,10,130/3,230/3,70,35]; % target positions
step_size = 1;
% manually pick the way pts for simulated human. change wap points directly 
% in getWayPts.m
inPara_gwp = struct('tar_pos',tar_pos,'scale',scale,'x_dis',x_dis,'type','h');
h_way_pts = getWayPts(inPara_gwp);

% apply different human speed between different way points.
% h_v = [2,3,1,1,1,1,1,1,1,1,3,1.5,2,3,2,1.5,4];
% apply different acceleration for human speed
% h_acl = -h.maxA+2*h.maxA*rand(1,300);
%%%

% define vertices of obstacles. Here use round and rectangular obstacles
% circle
c_set = [220/3-x_dis;40/3]; % coordinate of circle center
r_set = 10/3; % radius of circle
theta_set = {0:pi/8:2*pi};
% rectangle
rec_vert = {[0,36-x_dis;50,55];[44-x_dis,58-x_dis;25,30];[52-x_dis,68-x_dis;50,55];[68-x_dis,100-x_dis;30,34];[84-x_dis,100-x_dis;55,60]}; % vertices of rectangles.[low-left,upper-right]
inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta_set',{theta_set},...
    'rec_vert',{rec_vert},'type','obs');
obs_pos = getWayPts(inPara_gwp);

% moving obstacles
obs1.pos = [42;32];
obs1.v = 1;
obs1.s_time = 123;
obs1.dur = 30;
obs1.traj = obs1.pos;
obs2.pos = [53;82];
obs2.v = 1;
obs2.s_time = 195;
obs2.dur = 30;
obs2.traj = obs2.pos;

campus = field([xMin xMin+xLength/step_size yMin yMin+yLength/step_size],tar_pos);
campus.agentNum = 2;
% need to change this part and use ellipsoids to represent rectangle
% obstacles
ell_set = ellipseBound(rec_vert); %set of ellipses: [c;a;b]: c is the center of ellipsoid, a and b are the half length of long and short axis
campus.obs_info = {[c_set;r_set];ell_set}; % gives the center and radius of each round obstacle and vertices of rectangular obstacles

%% Simulation
% simulation parameters
kf = 350; % simulation length (/s)
agents = [h r];
hor = 5; % MPC horizon 
pre_type = 'IMM';%'extpol'; % 'extpol','IMM'. specify the method for predicting human motion
plan_type = 'greedy1'; % 'MPC','greedy1','greedy0'. specify the method for robot controller
samp_rate = 20; % sampling rate (/Hz)
safe_dis = 2; %safe distance between human and robot
safe_marg = 2; % safety margin between human the the obstacle
mpc_dt = 0.5; % sampling time for model discretization used in MPC

% initialize variables
obv_traj = zeros(3,0); % observed human trajectory; first row denotes the time [t,x,y]
obs1_obv_traj = zeros(3,0); % observed human trajectory; first row denotes the time [t,x,y]
obs2_obv_traj = zeros(3,0); % observed human trajectory; first row denotes the time [t,x,y]
est_state = zeros(4,mpc_dt*samp_rate,kf); % estimated human states for every second [x,vx,y,vy];
pre_traj = zeros(2,hor+1,kf); % current and predicted future human trajectory [x,y]
plan_state = zeros(4,hor+1,kf); % robot's current and planned future state [x,y,v]
r_state = [[r.currentPos;r.currentV],zeros(4,kf-1)]; % robot's actual state [x,y,theta,v]
r_input = zeros(2,kf); % robot's actual control input [w,a]
wp_cnt = 1; % the waypoint that the human is heading for
h_tar_wp = h_way_pts(:,wp_cnt); % the way point that the human is heading for

% the following snippet seems not in use, but only for purpose of the input
% of agentMove for robot
addpath('./sim_res')
load('x_pos_pre_imm','x_pos_pre_imm')
load('y_pos_pre_imm','y_pos_pre_imm')
pos_pre_imm = zeros(2,size(x_pos_pre_imm,1)-1,size(x_pos_pre_imm,2));
for ii = 1:size(x_pos_pre_imm,2)
    pos_pre_imm(:,:,ii) = [x_pos_pre_imm(2:end,ii)';y_pos_pre_imm(2:end,ii)'];
end
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
    parameter_IMM_UKF;  % parameter for IMM
    inPara_ams = struct('campus',campus,'agents',agents,'h_tar_wp',h_tar_wp,...
        'obv_traj',obv_traj1','est_state',est_state,...
        'pre_traj',pre_traj,'plan_state',plan_state,'r_state',r_state,'r_input',r_input,...
        'k',k,'hor',hor,'pre_type',pre_type,'samp_rate',samp_rate,...
        'safe_dis',safe_dis,'mpc_dt',mpc_dt,'safe_marg',safe_marg,...
        'agentIndex',agentIndex,'plan_type',plan_type,'samp_num',samp_num,...
        'pos_pre_imm',pos_pre_imm);
    [outPara_ams] = agentMove(inPara_ams);
    agents = outPara_ams.agents;
    est_state = outPara_ams.est_state;
    pre_traj = outPara_ams.pre_traj;
    plan_state = outPara_ams.plan_state; 
    r_state = outPara_ams.r_state;
    r_input = outPara_ams.r_input;
    %}
    
    %% obstacle moves
    % 1st obs moves horizontally
    if (k < obs1.s_time)
        cur_pos = obs1.pos;
        next_pos = obs1.pos;
    elseif (k >= obs1.s_time) && (k<= obs1.s_time+obs1.dur-1)
        cur_pos = obs1.pos;
        obs1.pos = obs1.pos+[0;obs1.v]*mpc_dt;
        next_pos = obs1.pos;
    else
        cur_pos = obs1.pos;
        next_pos = obs1.pos;
    end
    obs1.traj = [obs1.traj,obs1.pos];
    
    for ii = 1:samp_num
        % observed obs1 position
        tmp_t = (k-1)*mpc_dt+(ii-1)/samp_rate;
        obs1_obv_traj = [obs1_obv_traj,[tmp_t;cur_pos+(next_pos-cur_pos)*(ii-1)/samp_num]];
    end
    %}
    
    % 2nd obs moves vertically (along -y direction)
    if (k < obs2.s_time)
        cur_pos = obs2.pos;
        next_pos = obs2.pos;
    elseif (k >= obs2.s_time) && (k<= obs2.s_time+obs2.dur-1)
        cur_pos = obs2.pos;
        obs2.pos = obs2.pos-[0;obs2.v]*mpc_dt;
        next_pos = obs2.pos;
    else
        cur_pos = obs2.pos;
        next_pos = obs2.pos;
    end
    obs2.traj = [obs2.traj,obs2.pos];
    
    for ii = 1:samp_num
        % observed obs2 position
        tmp_t = (k-1)*mpc_dt+(ii-1)/samp_rate;
        obs2_obv_traj = [obs2_obv_traj,[tmp_t;cur_pos+(next_pos-cur_pos)*(ii-1)/samp_num]];
    end
    %% plot trajectories
    % Plot future agent positions
    
    % plot specifications
    color_agent = {'r','g','k','b','m','b'};
    marker_agent = {'o','^','*','d'};
    line_agent = {'-','-','-','-'};
    line_wid = {2,2,2,2};
    orange = [1 204/255 0];
    color_target = {'m','b',orange};
    figure;
    hold on

    % draw targets
    for jj = 1:campus.targetNum
        h = plot(campus.targetPos(1,jj),campus.targetPos(2,jj),'MarkerSize',15);
        set(h,'Marker','p');
        set(h,'linewidth',2);
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
        h1 = plot(tmp_agent.traj(1,:),tmp_agent.traj(2,:),'markers',2);
        set(h1,'MarkerFaceColor',color_agent{ii});
        set(h1,'MarkerEdgeColor',color_agent{ii});
        set(h1,'Color',color_agent{ii});
        set(h1,'LineStyle',line_agent{ii});
        set(h1,'Marker',marker_agent{ii});
        set(h1,'linewidth',line_wid{ii});
        h2 = plot(tmp_agent.currentPos(1),tmp_agent.currentPos(2),color_agent{ii},'markers',2);
        set(h2,'MarkerFaceColor',color_agent{ii});
        set(h2,'MarkerEdgeColor',color_agent{ii});
        set(h2,'Color',color_agent{ii});
        set(h2,'LineStyle',line_agent{ii});
        set(h2,'Marker',marker_agent{ii});
    end
    % predicted human positions
    h3 = plot(pre_traj(1,:,k),pre_traj(2,:,k),color_agent{3},'markers',2);
    set(h3,'MarkerFaceColor',color_agent{3});
    set(h3,'MarkerEdgeColor',color_agent{3});
    set(h3,'Color',color_agent{3});
    set(h3,'LineStyle',line_agent{3});
    set(h3,'Marker',marker_agent{3});
    set(h3,'linewidth',line_wid{3});
    hold on
    % I don't know why I need to repeat this snippet here.
    %{ 
    c_set = [pre_traj(1,2:end,k);pre_traj(2,2:end,k)];
    r_set = [safe_dis/2;safe_dis/2;safe_dis/2;safe_dis/2;safe_dis/2];
    theta_set1 = {{0:pi/8:2*pi;0:pi/8:2*pi;0:pi/8:2*pi;0:pi/8:2*pi;0:pi/8:2*pi}};
    inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta_set',theta_set1,'type','obs');
    circle_pos = getWayPts(inPara_gwp);
    for jj = 1:size(circle_pos)
        tmp_pos = circle_pos{jj};
        plot(tmp_pos(1,:),tmp_pos(2,:));%,color_agent{3});
    end
    %}
    
    % planned robot trajectory
    h4 = plot(plan_state(1,:,k),plan_state(2,:,k),color_agent{4},'markers',2);
    set(h4,'MarkerFaceColor',color_agent{4});
    set(h4,'MarkerEdgeColor',color_agent{4});
    set(h4,'Color',color_agent{4});
    set(h4,'LineStyle',line_agent{4});
    set(h4,'Marker',marker_agent{4});
    set(h1,'linewidth',line_wid{4});
    % I don't know why I need to repeat this snippet here.
    %{
    c_set = [plan_state(1,2:end,k);plan_state(2,2:end,k)];
    r_set = [safe_dis/2;safe_dis/2;safe_dis/2;safe_dis/2;safe_dis/2];
    theta_set1 = {{0:pi/8:2*pi;0:pi/8:2*pi;0:pi/8:2*pi;0:pi/8:2*pi;0:pi/8:2*pi}};
    inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta_set',theta_set1,'type','obs');
    circle_pos = getWayPts(inPara_gwp);
    for jj = 1:size(circle_pos)
        tmp_pos = circle_pos{jj};
        plot(tmp_pos(1,:),tmp_pos(2,:));%,color_agent{3});
    end
    %}
    
    % moving obstacles' trajectory
    %
    if k >= obs1.s_time
        h5 = plot(obs1.traj(1,:),obs1.traj(2,:),color_agent{5},'markers',2);
    end
    %}
    if k >= obs2.s_time
        h6 = plot(obs2.traj(1,:),obs2.traj(2,:),color_agent{6},'markers',2);
    end    
    
    % refine plot
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

% pre_traj = pos_pre_imm;
%% save simulation result
%{
% save data
% if the data is a decimal, replace the '.' with 'p'
str_safe_dis = strrep(num2str(safe_dis),'.','p');
str_safe_marg = strrep(num2str(safe_marg),'.','p');
str_h_v = strrep(num2str(agents(1).currentV),'.','p');

str_t = datestr(clock,'dd-mmm-yyyy_HHMMSS');
folder_path = ('./sim_res');
data_name = sprintf('sim_traj_%s_%s_%s_%s_%s_%s.mat',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name = fullfile (folder_path,data_name);
save(file_name,'obv_traj','est_state','pre_traj','plan_state','r_state','r_input');

% save obstacle trajectories
data_name = sprintf('sim_obs1_traj_%s_%s_%s_%s_%s_%s.mat',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name = fullfile (folder_path,data_name);
save(file_name,'obs1_obv_traj');

data_name = sprintf('sim_obs2_traj_%s_%s_%s_%s_%s_%s.mat',...
    pre_type,plan_type,str_safe_dis,str_safe_marg,str_h_v,str_t);
file_name = fullfile (folder_path,data_name);
save(file_name,'obs2_obv_traj');

% save plot
folder_path = ('./sim_res');
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