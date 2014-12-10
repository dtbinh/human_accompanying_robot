% process simulation result
%% imm vs extpol
% load data
%
clear
addpath('.\sim_res');
load('h_state.mat');
sim_len = size(h_state,2); % simulation time
h_traj = h_state(1:2,:); % h_state contains [x;y;heading]
% safe_dist1 = 2; % safe_dist for imm
% safe_dist2 = 1.5; % safe_dist for extpol
% h_v = 2;

% imm
load('sim_traj_IMM_MPC_08-Dec-2014_2_2','pre_traj','r_state');
imm_pre_traj = pre_traj(:,:,1:sim_len); % predicted human position
% imm_r_pos = r_state(1:2,1:sim_len);
% imm_r_v = r_state(3,1:sim_len);
% extrapolation
% load('sim_traj_extpol_MPC_08-Dec-2014_1p5_2','pre_traj','r_state');
% extpol_pre_traj-09-Dec-2014 saves the predicted trajectory using
% extrapolation. this is a tempory way.
load('extpol_pre_traj-08-Dec-2014','pre_traj');
extpol_pre_traj = pre_traj(:,:,1:sim_len);
% extpol_r_pos = r_state(1:2,1:sim_len);
% extpol_r_v = r_state(3,1:sim_len);

hor = size(pre_traj,2)-1; % planning horizon

%% compare the prediction
pred_err = zeros(2,sim_len-1); % prediction error. 1st row for imm and 2nd row for extpol
for ii = 1:sim_len-1
    if ii+hor <= sim_len-1
        dif_vec1 = imm_pre_traj(:,:,ii)-h_traj(:,ii:ii+hor);
        dif_vec2 = extpol_pre_traj(:,:,ii)-h_traj(:,ii:ii+hor);
        pred_err(1,ii) = sum(sqrt(sum(dif_vec1.^2,1)))/(hor+1);
        pred_err(2,ii) = sum(sqrt(sum(dif_vec2.^2,1)))/(hor+1);
    else
        tmp_len = sim_len-ii;
        dif_vec1 = imm_pre_traj(:,1:tmp_len,ii)-h_traj(:,ii:sim_len-1);
        dif_vec2 = extpol_pre_traj(:,1:tmp_len,ii)-h_traj(:,ii:sim_len-1);
        pred_err(1,ii) = sum(sqrt(sum(dif_vec1.^2,1)))/tmp_len;
        pred_err(2,ii) = sum(sqrt(sum(dif_vec2.^2,1)))/tmp_len;
    end
end
ave_pred_err = mean(pred_err,2);
max_pred_err = max(pred_err,[],2);
figure
hold on
plot(1:size(pred_err,2),pred_err(1,:),'r')
plot(1:size(pred_err,2),pred_err(2,:),'b')
legend('imm','extpol')
title('prediction difference')
%}
%% compare the motion planning
% code is correct, but using the wrong data. motion planning comparison
% should compare mpc and greedy, not imm vs extrapolation
%{
pos_dif = zeros(2,sim_len-1);
v_dif = zeros(2,sim_len-1);

for ii = 1:sim_len-1
    pos_dif(1,ii) = sqrt(norm(imm_r_pos(:,ii)-h_traj(:,ii),2))-safe_dist1;
    pos_dif(2,ii) = sqrt(norm(extpol_r_pos(:,ii)-h_traj(:,ii),2))-safe_dist2;
end
ave_pos_dif = mean(pos_dif,2);
max_pos_dif = max(abs(pos_dif),[],2);
figure
hold on
plot(1:size(pos_dif,2),pos_dif(1,:),'r')
plot(1:size(pos_dif,2),pos_dif(2,:),'b')
legend('imm','extpol')
title('distance difference')

v_dif(1,:) = imm_r_v(1:end-1)-h_v;
v_dif(2,:) = extpol_r_v(1:end-1)-h_v;
ave_v_dif = mean(v_dif,2);
max_v_dif = max(v_dif,[],2);
figure
hold on
plot(1:size(v_dif,2),v_dif(1,:),'r')
plot(1:size(v_dif,2),v_dif(2,:),'b')
legend('imm','extpol')
title('velocity difference')
legend('imm','extpol')
%}
%% mpc vs greedy

%% compare the motion planning
% load data
addpath('.\sim_res');
load('h_state.mat');
sim_len = size(h_state,2); % simulation time
h_traj = h_state(1:2,:); % h_state contains [x;y;heading]
safe_dist = 2; % safe_dist for imm
h_v = 2;

% mpc
load('sim_traj_IMM_MPC_08-Dec-2014_2_2','pre_traj','r_state');
imm_pre_traj = pre_traj(:,:,1:sim_len); % predicted human position
imm_r_pos = r_state(1:2,1:sim_len);
imm_r_v = r_state(3,1:sim_len);
% greedy
load('sim_traj_IMM_greedy1_09-Dec-2014_2_2','pre_traj','r_state');
extpol_pre_traj = pre_traj(:,:,1:sim_len);
extpol_r_pos = r_state(1:2,1:sim_len);
extpol_r_v = r_state(3,1:sim_len);

pos_dif = zeros(2,sim_len-1);
v_dif = zeros(2,sim_len-1);

for ii = 1:sim_len-1
    dif_vec1 = imm_r_pos(:,ii)-h_traj(:,ii+1);
    dif_vec2 = extpol_r_pos(:,ii)-h_traj(:,ii+1);
    pos_dif(1,ii) = norm(dif_vec1,2)-safe_dist;
    pos_dif(2,ii) = norm(dif_vec2,2)-safe_dist;
end
ave_pos_dif = mean(pos_dif,2);
max_pos_dif = max(abs(pos_dif),[],2);
figure
hold on
plot(1:size(pos_dif,2),pos_dif(1,:),'r')
plot(1:size(pos_dif,2),pos_dif(2,:),'b')
legend('mpc','greedy')
title('distance difference')

v_dif(1,:) = imm_r_v(1:end-1)-h_v;
v_dif(2,:) = extpol_r_v(1:end-1)-h_v;
ave_v_dif = mean(v_dif,2);
max_v_dif = max(abs(v_dif),[],2);
figure
hold on
plot(1:size(v_dif,2),v_dif(1,:),'r')
plot(1:size(v_dif,2),v_dif(2,:),'b')
legend('mpc','greedy')
title('velocity difference')