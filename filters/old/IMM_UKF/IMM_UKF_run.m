function [x_est,y_est,x_pos_pre,y_pos_pre] = IMM_UKF_run()

% Run Simulink Part
sim('IMM_UKF.slx');
end