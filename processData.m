% process simulation result
%% imm vs single-model
% load data
%{
clear
addpath('./sim_res');
load('h_state_1p5.mat');
sim_len = size(h_state,2); % simulation time
h_traj = h_state(1:2,:); % h_state contains [x;y;heading]
% safe_dist1 = 2; % safe_dist for imm
% safe_dist2 = 1.5; % safe_dist for extpol
% h_v = 2;

% % extrapolation
% % load('sim_traj_extpol_MPC_08-Dec-2014_1p5_2','pre_traj','r_state');
% % extpol_pre_traj-09-Dec-2014 saves the predicted trajectory using
% % extrapolation. this is a tempory way.
% load('sim_traj_extpol_greedy1_2_2_1p5_09-Dec-2014_204347','pre_traj');
% extpol_pre_traj = pre_traj(:,:,1:sim_len);
% % extpol_r_pos = r_state(1:2,1:sim_len);
% % extpol_r_v = r_state(3,1:sim_len);

% single model
load('x_pos_LKF','x_pos_LKF')
load('y_pos_LKF','y_pos_LKF')

lkf_pre_traj = zeros(2,size(x_pos_LKF,1)-2,size(x_pos_LKF,2));
lkf_pre_traj(1,:,:) = x_pos_LKF(3:end,:);
lkf_pre_traj(2,:,:) = y_pos_LKF(3:end,:);

load('x_pos_UKF','x_pos_UKF')
load('y_pos_UKF','y_pos_UKF')

ukf_pre_traj = zeros(2,size(x_pos_UKF,1)-2,size(x_pos_UKF,2));
ukf_pre_traj(1,:,:) = x_pos_UKF(3:end,:);
ukf_pre_traj(2,:,:) = y_pos_UKF(3:end,:);

% imm
load('x_pos_imm_LKF','x_pos_imm_LKF')
load('y_pos_imm_LKF','y_pos_imm_LKF')

imm_lkf_pre_traj = zeros(2,size(x_pos_imm_LKF,1)-2,size(x_pos_imm_LKF,2));
imm_lkf_pre_traj(1,:,:) = x_pos_imm_LKF(3:end,:);
imm_lkf_pre_traj(2,:,:) = y_pos_imm_LKF(3:end,:);

load('x_pos_imm_UKF','x_pos_imm_UKF')
load('y_pos_imm_UKF','y_pos_imm_UKF')

imm_ukf_pre_traj = zeros(2,size(x_pos_imm_UKF,1)-2,size(x_pos_imm_UKF,2));
imm_ukf_pre_traj(1,:,:) = x_pos_imm_UKF(3:end,:);
imm_ukf_pre_traj(2,:,:) = y_pos_imm_UKF(3:end,:);

hor = size(imm_lkf_pre_traj,2);

% compare the prediction
%
pred_err = zeros(4,sim_len-1); % prediction error. 1st row: lfk; 2nd: ukf; 3rd: imm_lkf; 4th: imm_ukf
for ii = 1:sim_len-2
    if ii+hor <= sim_len-1
        dif_vec1 = lkf_pre_traj(:,:,ii)-h_traj(:,ii+1:ii+hor);
        dif_vec2 = ukf_pre_traj(:,:,ii)-h_traj(:,ii+1:ii+hor);
        dif_vec3 = imm_lkf_pre_traj(:,:,ii)-h_traj(:,ii+1:ii+hor);
        dif_vec4 = imm_ukf_pre_traj(:,:,ii)-h_traj(:,ii+1:ii+hor);
        pred_err(1,ii) = sum(sqrt(sum(dif_vec1.^2,1)))/(hor);
        pred_err(2,ii) = sum(sqrt(sum(dif_vec2.^2,1)))/(hor);
        pred_err(3,ii) = sum(sqrt(sum(dif_vec3.^2,1)))/(hor);
        pred_err(4,ii) = sum(sqrt(sum(dif_vec4.^2,1)))/(hor);
    else
        tmp_len = sim_len-ii-1;
        dif_vec1 = lkf_pre_traj(:,1:tmp_len,ii)-h_traj(:,ii+1:sim_len-1);
        dif_vec2 = ukf_pre_traj(:,1:tmp_len,ii)-h_traj(:,ii+1:sim_len-1);
        dif_vec3 = imm_lkf_pre_traj(:,1:tmp_len,ii)-h_traj(:,ii+1:sim_len-1);
        dif_vec4 = imm_ukf_pre_traj(:,1:tmp_len,ii)-h_traj(:,ii+1:sim_len-1);
        pred_err(1,ii) = sum(sqrt(sum(dif_vec1.^2,1)))/(tmp_len);
        pred_err(2,ii) = sum(sqrt(sum(dif_vec2.^2,1)))/(tmp_len);
        pred_err(3,ii) = sum(sqrt(sum(dif_vec3.^2,1)))/(tmp_len);
        pred_err(4,ii) = sum(sqrt(sum(dif_vec4.^2,1)))/(tmp_len);
    end
end
ave_pred_err = mean(pred_err,2);
std_pred_err = (std(pred_err'))';
max_pred_err = max(pred_err,[],2);
h1 = figure;
box on
hold on
plot((1:size(pred_err,2))*0.5,pred_err(1,:),'m--','LineWidth',2)
plot((1:size(pred_err,2))*0.5,pred_err(3,:),'r:','LineWidth',2)
plot((1:size(pred_err,2))*0.5,pred_err(2,:),'c--','LineWidth',2)
plot((1:size(pred_err,2))*0.5,pred_err(4,:),'b-.','LineWidth',2)
legend('lkf','imm-lkf','ukf','imm-ukf')
grid on
title('Prediction difference')
xlabel('time [sec]')
ylabel('distance [m]')
xlim([0,160])
saveas(h1,'imm_vs_single','fig')
% fig2Pdf('imm_vs_single',300,h1)
saveTightFigure(h1,'imm_vs_single')
%}

%% mpc vs reactive
%{
% compare the motion planning
% load data
%
str_t = datestr(clock,'dd-mmm-yyyy_HHMMSS');
addpath('./sim_res');
load('human_traj.mat'); %h_state_1p5
sim_len = size(human_traj,2); % simulation time
h_traj = human_traj(1:2,:); % h_state contains [x;y;heading]
safe_dist = 1; % safe_dist for imm
h_v = 1.5;

% mpc
load('sim_traj_UKF_MPC_1_1_1p5_04-Jun-2015_200103','pre_traj','r_state','plan_state');
mpc_pre_traj = pre_traj(:,:,1:sim_len); % predicted human position
mpc_r_plan_pos = plan_state(1:2,:,:);
mpc_r_pos = r_state(1:2,1:sim_len);
mpc_r_v = r_state(4,1:sim_len);
% greedy
%{
load('sim_traj_IMM_greedy1_2_2_1p5_25-Apr-2015_215112','r_state','plan_state');
grd_pre_traj = pre_traj(:,:,1:sim_len);
grd_r_plan_pos = plan_state(1:2,:,:);
grd_r_pos = r_state(1:2,1:sim_len);
grd_r_v = r_state(4,1:sim_len);
%}

 % position difference between human and robot. 
 % Four elements: 
 % 1. diff between actual positions of h and r for mpc
 % 2. diff between actual positions of h and r for greedy
 % 3. diff between predicted positions of h and actual positions of r for mpc
 % 4. diff between predicted positions of h and actual positions of r for mpc
pos_dif = zeros(4,sim_len-1);
% velocity difference
v_dif = zeros(2,sim_len-1);
mpc_plan_pos_dif = zeros(6,sim_len-1); % difference between human estimated positon and robot planned positions in the prediction horizon
% grd_plan_pos_dif = zeros(6,sim_len-1); % difference between human estimated positon and robot planned positions in the prediction horizon

for ii = 1:sim_len-2
    dif_vec1 = mpc_r_pos(:,ii)-h_traj(:,ii);
%     dif_vec2 = grd_r_pos(:,ii)-h_traj(:,ii);
    pos_dif(1,ii) = norm(dif_vec1,2)-safe_dist;
%     pos_dif(2,ii) = norm(dif_vec2,2)-safe_dist;
    dif_vec3 = mpc_r_pos(:,ii)-mpc_pre_traj(:,1,ii);
%     dif_vec4 = grd_r_pos(:,ii)-grd_pre_traj(:,1,ii);
    pos_dif(3,ii) = norm(dif_vec3,2)-safe_dist;
%     pos_dif(4,ii) = norm(dif_vec4,2)-safe_dist;
    dif_vec5 = mpc_r_plan_pos(:,:,ii)-mpc_pre_traj(:,:,ii);
%     dif_vec6 = grd_r_plan_pos(:,2:end,ii)-grd_pre_traj(:,2:end,ii);
    mpc_plan_pos_dif(:,ii) = sqrt((sum(dif_vec5.*dif_vec5,1))');
%     grd_plan_pos_dif(:,ii) = sqrt((sum(dif_vec6.*dif_vec6,1))');
end


ave_pos_dif = mean(pos_dif,2);
max_pos_dif = max(abs(pos_dif),[],2);
h2 = figure;
hold on
box on
% plot((1:size(pos_dif,2))*0.5,pos_dif(2,:),'r','LineWidth',2)
plot((1:size(pos_dif,2))*0.5,pos_dif(1,:),'b','LineWidth',2)
% legend('reactive','mpc')
grid on
title('Distance difference')
xlabel('time [sec]')
ylabel('distance [m]')
xlim([0,160])
% saveas(h2,'mpc_vs_reac_dis_diff','fig')
% fig2Pdf('mpc_vs_reac_dis_diff',300,h2)
saveas(h2,sprintf('sim_res/ukf_mpc_dis_diff_%s',str_t),'fig')
fig2Pdf(sprintf('sim_res/ukf_mpc_dis_diff_%s',str_t),300,h2)


%{
h3 = figure;
hold on
box on
plot((1:size(pos_dif,2))*0.5,pos_dif(4,:),'r--','LineWidth',2)
plot((1:size(pos_dif,2))*0.5,pos_dif(3,:),'b','LineWidth',2)
legend('reactive','mpc')
grid on
title('Distance difference between estimated human positions and robot positions')
xlabel('time [sec]')
ylabel('distance [m]')
xlim([0,160])
saveas(h3,'mpc_vs_reac_dis_diff','fig')
%}

v_dif(1,:) = mpc_r_v(1:end-1)-h_v;
% v_dif(2,:) = grd_r_v(1:end-1)-h_v;
ave_v_dif = mean(v_dif,2);
max_v_dif = max(abs(v_dif),[],2);
h4 = figure;
hold on
box on
% plot((1:size(v_dif,2))*0.5,v_dif(2,:),'r','LineWidth',2)
plot((1:size(v_dif,2))*0.5,v_dif(1,:),'b','LineWidth',2)
% legend('reactive','mpc')
grid on
title('Velocity difference')
xlabel('time[sec]')
ylabel('speed [m/s]')
xlim([0,160])
% saveas(h4,'mpc_vs_reac_vel_diff','fig')
% fig2Pdf('mpc_vs_reac_vel_diff',300,h4)
saveas(h4,sprintf('sim_res/ukf_mpc_vel_diff_%s',str_t),'fig')
fig2Pdf(sprintf('sim_res/ukf_mpc_vel_diff_%s',str_t),300,h4)
 
% save('sim_res.mat','ave_pred_err','max_pred_err',...
%     'ave_pos_dif','max_pos_dif','ave_v_dif','max_v_dif');
%}

%% IMM-UKF vs UKF vs no-prediction (mpc as the planner)
%
% compare different prediction method using the same planner
% load data
%
str_t = datestr(clock,'dd-mmm-yyyy_HHMMSS');
addpath('./sim_res');
load('human_traj.mat'); %h_state_1p5
sim_len = size(human_traj,2); % simulation time
h_traj = human_traj(1:2,:); % h_state contains [x;y;heading]
safe_dist = 1; % safe_dist for imm
h_v = 1.5;

% imm-ukf
load('sim_traj_IMM-UKF_MPC_1_1_1p5_21-Jun-2015_155124','pre_traj','r_state','plan_state');
imm_pre_traj = pre_traj(:,:,1:sim_len); % predicted human position
imm_r_plan_pos = plan_state(1:2,:,:);
imm_r_pos = r_state(1:2,1:sim_len);
imm_r_v = r_state(4,1:sim_len);
% ukf
%
load('sim_traj_UKF_MPC_1_1_1p5_21-Jun-2015_152805','pre_traj','r_state','plan_state');
ufk_pre_traj = pre_traj(:,:,1:sim_len);
ufk_r_plan_pos = plan_state(1:2,:,:);
ufk_r_pos = r_state(1:2,1:sim_len);
ufk_r_v = r_state(4,1:sim_len);
%}
% no-prediction
load('sim_traj_IMM-UKF_greedy0_1_1_1p5_21-Jun-2015_163734','pre_traj','r_state','plan_state');
nop_pre_traj = pre_traj(:,:,1:sim_len);
nop_r_plan_pos = plan_state(1:2,:,:);
nop_r_pos = r_state(1:2,1:sim_len);
nop_r_v = r_state(4,1:sim_len);

 % position difference between human and robot. 
 % Six elements: 
 % 1. diff between actual positions of h and r for imm-ukf
 % 2. diff between actual positions of h and r for ukf
 % 3. diff between actual positions of h and r for no-prediction
 % 4. diff between predicted positions of h and actual positions of r for imm-ukf
 % 5. diff between predicted positions of h and actual positions of r for ukf
 % 6. diff between predicted positions of h and actual positions of r for no-prediction
pos_dif = zeros(6,sim_len-1);
% velocity difference
v_dif = zeros(3,sim_len-1);
imm_plan_pos_dif = zeros(6,sim_len-1); % difference between human estimated positon and robot planned positions in the prediction horizon
ukf_plan_pos_dif = zeros(6,sim_len-1); % difference between human estimated positon and robot planned positions in the prediction horizon
nop_plan_pos_dif = zeros(6,sim_len-1); % difference between human estimated positon and robot planned positions in the prediction horizon

for ii = 1:sim_len-2
    dif_vec1 = imm_r_pos(:,ii)-h_traj(:,ii);
    dif_vec2 = ufk_r_pos(:,ii)-h_traj(:,ii);
    dif_vec3 = nop_r_pos(:,ii)-h_traj(:,ii);
    pos_dif(1,ii) = norm(dif_vec1,2);
    pos_dif(2,ii) = norm(dif_vec2,2);
    pos_dif(3,ii) = norm(dif_vec3,2);
    dif_vec4 = imm_r_pos(:,ii)-imm_pre_traj(:,1,ii);
    dif_vec5 = ufk_r_pos(:,ii)-ufk_pre_traj(:,1,ii);
    dif_vec6 = nop_r_pos(:,ii)-nop_pre_traj(:,1,ii);
    pos_dif(4,ii) = norm(dif_vec4,2);
    pos_dif(5,ii) = norm(dif_vec5,2);
    pos_dif(6,ii) = norm(dif_vec6,2);
    dif_vec7 = imm_r_plan_pos(:,:,ii)-imm_pre_traj(:,:,ii);
    dif_vec8 = ufk_r_plan_pos(:,:,ii)-ufk_pre_traj(:,:,ii);
    dif_vec9 = nop_r_plan_pos(:,:,ii)-nop_pre_traj(:,:,ii);
    imm_plan_pos_dif(:,ii) = sqrt((sum(dif_vec7.*dif_vec7,1))');
    ukf_plan_pos_dif(:,ii) = sqrt((sum(dif_vec8.*dif_vec8,1))');
    nop_plan_pos_dif(:,ii) = sqrt((sum(dif_vec9.*dif_vec9,1))');
end

ave_pos_dif = mean(pos_dif,2);
std_pos_dif = std(pos_dif,0,2);
max_pos_dif = max(abs(pos_dif),[],2);
h2 = figure;
hold on
box on
plot((1:size(pos_dif,2))*0.5,pos_dif(1,:),'b','LineWidth',2)
plot((1:size(pos_dif,2))*0.5,pos_dif(2,:),'r','LineWidth',2)
plot((1:size(pos_dif,2))*0.5,pos_dif(3,:),'k','LineWidth',2)
legend('IMM-UKF','UKF','NOP')
grid on
title('Distance difference')
xlabel('time [sec]')
ylabel('distance [m]')
xlim([0,160])
saveas(h2,sprintf('sim_res/imm_vs_ukf_vs_no-pred_dis_diff_%s',str_t),'fig')
fig2Pdf(sprintf('sim_res/imm_vs_ukf_vs_no-pred_dis_diff_%s',str_t),300,h2)
% saveas(h2,sprintf('sim_res/ukf_mpc_dis_diff_%s',str_t),'fig')
% fig2Pdf(sprintf('sim_res/ukf_mpc_dis_diff_%s',str_t),300,h2)


%{
h3 = figure;
hold on
box on
plot((1:size(pos_dif,2))*0.5,pos_dif(4,:),'r--','LineWidth',2)
plot((1:size(pos_dif,2))*0.5,pos_dif(3,:),'b','LineWidth',2)
legend('reactive','mpc')
grid on
title('Distance difference between estimated human positions and robot positions')
xlabel('time [sec]')
ylabel('distance [m]')
xlim([0,160])
saveas(h3,'mpc_vs_reac_dis_diff','fig')
%}

v_dif(1,:) = imm_r_v(1:end-1)-h_v;
v_dif(2,:) = ufk_r_v(1:end-1)-h_v;
v_dif(3,:) = nop_r_v(1:end-1)-h_v;
ave_v_dif = mean(v_dif,2);
std_v_dif = std(v_dif,0,2);
max_v_dif = max(abs(v_dif),[],2);
h4 = figure;
hold on
box on
plot((1:size(v_dif,2))*0.5,v_dif(1,:),'b','LineWidth',2)
plot((1:size(v_dif,2))*0.5,v_dif(2,:),'r','LineWidth',2)
plot((1:size(v_dif,2))*0.5,v_dif(3,:),'k','LineWidth',2)
legend('IMM-UKF','UKF','NOP')
grid on
title('Velocity difference')
xlabel('time[sec]')
ylabel('speed [m/s]')
xlim([0,160])
saveas(h4,sprintf('sim_res/imm_vs_ukf_vs_no-pred_vel_diff_%s',str_t),'fig')
fig2Pdf(sprintf('sim_res/imm_vs_ukf_vs_no-pred_vel_diff_%s',str_t),300,h4)
% saveas(h4,sprintf('sim_res/ukf_mpc_vel_diff_%s',str_t),'fig')
% fig2Pdf(sprintf('sim_res/ukf_mpc_vel_diff_%s',str_t),300,h4)
 
% save('sim_res.mat','ave_pred_err','max_pred_err',...
%     'ave_pos_dif','max_pos_dif','ave_v_dif','max_v_dif');
%}