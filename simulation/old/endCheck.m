% 11/24/14
function [outPara] = endCheck(inPara)
% Game over if human has gone to all the way points
% get input arguments
    h = inPara.h;    
    way_pts = inPara.h_way_pts;
    wp_cnt = inPara.wp_cnt;
    dt = inPara.mpc_dt;
    
    v = h.currentV;
    h_tar_wp = way_pts(:,wp_cnt);
    
    game_end = 0; 
    arr_wp_flag = 0;
    
    if norm (h.currentPos(1:2)-h_tar_wp,2) <= v*dt
        if wp_cnt == size(way_pts,2) % if the human has arrived at the last way point
            sprintf('simulation completes!')
            game_end = 1;
            arr_wp_flag = 1;
        else
            sprintf('human has arrived at %d th waypoint',wp_cnt)
            arr_wp_flag = 1;
        end
    end
      
    outPara = struct('game_end',game_end,'arr_wp_flag',arr_wp_flag);    
end         