% Parameter for IMM-LKF

T = 0.05; % for estimation [s]
T1 = 0.5; % for Prediction [s]
noise = 0; % Noise Variance [m]
SimulStop = length(obv_traj1)*T;


%% Plant

r2 = 0.1; % for Traj 3 (Not rapid repsonse but noisy)
% r2 = 0.1;

A1 = [1 T 0 0;...
      0 1 0 0;...
      0 0 1 T;...
      0 0 0 1];
 
A2 = [1 sin(r2*T)/r2 0 -(1-cos(r2*T))/r2;...
     0 cos(r2*T) 0 -sin(r2*T);...
     0 (1-cos(r2*T))/r2 1 sin(r2*T)/r2;...
     0 sin(r2*T) 0 cos(r2*T)];
 
A3 = [1 T1 0 0;...
      0 1 0 0;...
      0 0 1 T1;...
      0 0 0 1];

A4 = [1 sin(r2*T1)/r2 0 -(1-cos(r2*T1))/r2;...
     0 cos(r2*T1) 0 -sin(r2*T1);...
     0 (1-cos(r2*T1))/r2 1 sin(r2*T1)/r2;...
     0 sin(r2*T1) 0 cos(r2*T1)];
 
Bw = [T 0;1 0;0 T;0 1];
%Bw = [T^2/2 0;T 0;0 T^2/2;0 T];

C = [1 0 0 0; 0 0 1 0];
Dw = eye(2);

%% KF 
meas_noise = 1.5; 
process_noise = 1.5*10^(-2);  % for Traj 3 (NOT rapid repsonse but noisy)

%%% Original
Sv1 = meas_noise^2;
Sv2 = meas_noise^2;
Sw1 = process_noise^2;
Sw2 = process_noise^2;

% Sv1 = meas_noise1^2;
% Sw1 = process_noise^2;
% Sw2 = process_noise^2;

W1 = Bw*Sw1*Bw';
W2 = Bw*Sw2*Bw';
V1 = Dw*Sv1*Dw';
V2 = Dw*Sv2*Dw';


%Markov chain transition matrix

PHI=[0.99 0.01; 0.01 0.99];  %%%original
%PHI=[0.97 0.03; 0.03 0.97];


%Initial Condition
xhat_init=[obv_traj1(1,2);10^(-5) ;obv_traj1(1,3);0]; 
%xhat_init=[0;0 ;0;0]; 

P_init=eye(4)*1;
mu1_init=0.5;
mu2_init=0.5;


