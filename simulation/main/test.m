% scripts for debugging and testing

%% test the IMM
%{
clc; clear; close all;
load('obv_traj3_w_time.mat')% Load Tracjectory of Human
obv_traj1=obv_traj1';
parameter;  %%% Parameter File upload
[x_est,y_est,x_pos_pre,y_pos_pre] = IMM_Com_run(); %%% Run Simulink
% !IMM_2model_dhl_verN_ConvertTest;
% load IMM_2model_dhl_verN_ConvertTest.MAT
%}

%% test ROS topics
sub_r = rossubscriber('/robotState','std_msgs/Float64MultiArray');
r_msg = receive(sub_r,10)'
display(r_msg.Data)
