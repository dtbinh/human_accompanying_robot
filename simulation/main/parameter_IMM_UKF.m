% Parameter for IMM-UKF
% by D.H.Lee

T = 0.05; % for estimation [s]
T1 = 0.5; % for Prediction [s]
noise = 0.0; % Noise Variance [m]
SimulStop = length(obv_traj1)*T;


% Markov chain transition matrix
PHI=[0.99 0.01; 0.01 0.99];  %%%original
%PHI=[0.97 0.03; 0.03 0.97];


%Initial Condition
xhat_init = [obv_traj1(1,2);10^(-5);obv_traj1(1,3);0;10^(-5)];
%P_init=eye(5)*1;
mu1_init=0.5;
mu2_init=0.5;


