% addpath('.\IMM_module');
load('obv_traj4_w_time.mat')% Load Tracjectory of Human
obv_traj1=obv_traj1';
% parameter;  %%% Parameter File upload
outPara_gsp = getSimPara(obv_traj1);  %%% Parameter File upload
% SimulStop = outPara_gsp.SimulStop;
% A = outPara_gsp.A;
% A1 = outPara_gsp.A1;
% C = outPara_gsp.C;
% W1 = outPara_gsp.W1;
% W2 = outPara_gsp.W2;
% V1 = outPara_gsp.V1;
% PHI = outPara_gsp.PHI;
% xhat_init = outPara_gsp.xhat_init;
% P_init = outPara_gsp.P_init;
% mu1_init = outPara_gsp.mu1_init;
% mu2_init = outPara_gsp.mu2_init;
% est_state([1,2],k),est_state([3,4],k),pre_traj(1,:,k),pre_traj(2,:,k)
[x_est,y_est,x_pos_pre,y_pos_pre]= IMM_Com_run(outPara_gsp); %%% Run Simulink