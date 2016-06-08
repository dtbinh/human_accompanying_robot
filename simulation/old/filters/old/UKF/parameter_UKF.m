% parameter for UKF
% 2015-01-17 by D.H.Lee

T = 0.05; % for estimation
T1 = 0.5; % for Prediction
noise = 0.0;
SimulStop = length(obv_traj1)*T;
%SimulStop = 50;
