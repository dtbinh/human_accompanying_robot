% 11/24/2014
% this file is for the ME 290J course project
clc 
clearvars
clearvars -global % clear global variables
close all

%% Setup
set(0,'DefaultFigureWindowStyle','docked');
%%% define agents %%%
% Human agent 1
h = agent('human');
h.currentPos = [30;10;0];%[290;30;0];%[3.5;15.5;pi/2]; % [x y heading]
h.maxV = 5;
h.dPsi = pi/2;
h.maxDPsi = pi;

% Robot agent 1
r = agent('robot');
r.currentPos = [30;20;0];%[310;30;0]; %[23.5;0.5;0]; % [x y heading]
r.maxV = 10;
r.dPsi = pi/2;
r.maxDPsi = pi/4;

%%% Set field %%%
xLength = 300; 
yLength = 300; 
xMin = 0;
yMin = 0;
tar_pos = [90,275,270,250,20;30,50,130,270,210]; % target positions
step_size = 1;
% manually pick the way pts for simulated human
way_pts = [tar_pos(:,1:2),[230,230;90,110],tar_pos(:,3:4),[210,190;230,230],tar_pos(:,5)];
campus = field([xMin xMin+xLength/step_size yMin yMin+yLength/step_size],tar_pos);
campus.agentNum = 2;

obs_pos = zeros(2,4,3); % vetices of obstacles
obs_pos(:,:,1) = [250,300,300,250;90,90,110,110];
obs_pos(:,:,2) = [210,210,190,190;250,300,300,250];
obs_pos(:,:,3) = [0,50,50,0;140,140,160,160];

%% Simulation
% simulation parameters
kf = 500;
agents = [h r];
hor = 5; % MPC horizon
pre_type = 'extpol'; % 'IMM' specify the method for predicting human motion

% initialize variables
obv_traj = zeros(2,kf); % observed human position
est_state = zeros(5,kf);
pre_traj = zeros(2,hor,kf);
plan_path = zeros(2,hor,kf);
r_path = zeros(2,kf);
wp_cnt = 1; % counter for waypoints
h_tar_wp = way_pts(:,wp_cnt); % the way point that the human is heading for


for k = 1:kf
    %% waypoint check
    % check if the human needs to choose the next way point
    display(k)
    inPara_ec = struct('h',agents(1),'way_pts',way_pts,'wp_cnt',wp_cnt);
    [outPara_ec] = endCheck(inPara_ec);

    game_end = outPara_ec.game_end;
    arr_wp_flag = outPara_ec.arr_wp_flag; % the human has arrived at a way point
    
    if game_end == 1
        sprintf('total search time is %d',k-1)
        break
    end
    
    if arr_wp_flag == 1
        wp_cnt = wp_cnt+1;
        h_tar_wp = way_pts(:,wp_cnt);
    end
    
    %% both agents move
    % human moves to a way point and robot predicts and plans its own path
    inPara_ams = struct('agents',agents,'h_tar_wp',h_tar_wp,...
        'obv_traj',obv_traj,'est_state',est_state,...
        'pre_traj',pre_traj,'plan_path',plan_path,'r_path',r_path,'k',k,...
        'hor',hor,'pre_type',pre_type);
    [outPara_ams] = agentMove(campus,inPara_ams);
    agents = outPara_ams.agents;
    obv_traj = outPara_ams.obv_traj; % observed human trajectory
    est_state = outPara_ams.est_state; % estimated human states
    pre_traj = outPara_ams.pre_traj; % predicted human trajectory
    plan_path = outPara_ams.plan_path; % robot's planned path
    r_path = outPara_ams.r_path; % robot's actual path
        
    %% plot trajectories
    % Plot future agent positions
    
    % plot specifications
    color_agent = {'r'; 'g'};
    marker_agent = {'o','^'};
    line_agent = {'-','-'};
    orange = [1 204/255 0];
    color_target = {'m','b',orange};
    figure;
    hold on

    % draw targets
    for jj = 1:campus.targetNum
        h = plot(campus.targetPos(1,jj),campus.targetPos(2,jj),'MarkerSize',15);
        set(h,'Marker','p');
    end
    
    % draw obstacles
    for jj = 1:size(obs_pos,3)
        fill(obs_pos(1,:,jj),obs_pos(2,:,jj),'b');
    end
    
    xlim([0,campus.endpoints(2)]);
    ylim([0,campus.endpoints(4)]);
    
    % draw agent trajectory
    for ii = 1%1:length(agents)
        tmp_agent = agents(ii);
        h1 = plot(tmp_agent.traj(1,:),tmp_agent.traj(2,:),'markers',5);
        set(h1,'MarkerFaceColor',color_agent{ii});
        set(h1,'MarkerEdgeColor',color_agent{ii});
        set(h1,'LineStyle',line_agent{ii});
        set(h1,'Marker',marker_agent{ii});
        h2 = plot(tmp_agent.currentPos(1),tmp_agent.currentPos(2),color_agent{ii},'markers',10);
        set(h2,'MarkerFaceColor',color_agent{ii});
        set(h2,'MarkerEdgeColor',color_agent{ii});
        set(h2,'LineStyle',line_agent{ii});
        set(h2,'Marker',marker_agent{ii});
    end
    grid minor
    set(gca,'MinorGridLineStyle','-','XColor',[0.5 0.5 0.5],'YColor',[0.5 0.5 0.5])
    axis equal
    xlim([0,xLength]);ylim([0,yLength]);
    %}
end
