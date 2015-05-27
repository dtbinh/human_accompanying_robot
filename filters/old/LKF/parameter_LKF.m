% parameter for LKF

T = 0.05; % for estimation
T1 = 0.5; % for Prediction
noise = 0.0;
SimulStop = length(obv_traj1)*T;

%% Plant

r1 = 0.1; % for Traj 3 (Not rapid repsonse but noisy)


 
A = [1 T 0 0;...
      0 1 0 0;...
      0 0 1 T;...
      0 0 0 1];
  

A1 = [1 T1 0 0;...
      0 1 0 0;...
      0 0 1 T1;...
      0 0 0 1];

Bw = [T 0;1 0;0 T;0 1];

C = [1 0 0 0; 0 0 1 0];
Dw = eye(2);

%% KF 
meas_noise1 = 1.5; 
process_noise = 1.5*10^(-2);  % for Traj 3 (NOT rapid repsonse but noisy)


Sv1 = meas_noise1^2;
Sw1 = process_noise^2;

W1 = Bw*Sw1*Bw';
V1 = Dw*Sv1*Dw';

% Initial Condition
xhat_init=[obv_traj1(1,2);10^(-5) ;obv_traj1(1,3);0]; 
%xhat_init=[obv_traj1(1,2);10 ;obv_traj1(1,3);0]; 

P_init=eye(4)*1;
mu1_init=0.5;
mu2_init=0.5;
