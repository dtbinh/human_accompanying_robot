% clc; clear;
load('obv_traj4_w_time.mat')% Load Tracjectory of Human
obv_traj1=obv_traj1';
parameter;  %%% Parameter File upload
[x_est,y_est,x_pos_pre,y_pos_pre] = IMM_Com_run(); %%% Run Simulink

