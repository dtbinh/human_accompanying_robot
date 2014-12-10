function [x_est,y_est,x_pos_pre,y_pos_pre,x_pos_est,input,time] = IMM_Com_run()
% 290J TermPJT IMM part
% by D.H.Lee
% Run Simulink Part
sim('IMM_2model_turn_rate_Com.slx');
end