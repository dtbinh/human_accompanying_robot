function way_pts = getWayPts(inPara)
type = inPara.type;
if (strcmp(type,'h'))
    tar_pos = inPara.tar_pos;
    x_dis = inPara.x_dis;
    % human way points
    % generate way points for human trajectory around the round obstacle
    center = [220/3-x_dis;40/3]; %center of circle
    r = 5; %radius
    theta = -3/4*pi:pi/8:pi/4;
    cor = zeros(2,length(theta)); %coordinates
    for ii = 1:length(theta)
        cor(1,ii) = center(1)+r*cos(theta(ii));
        cor(2,ii) = center(2)+r*sin(theta(ii));
    end
%     figure
%     hold on
%     fill(cor(1,:),cor(2,:),'r');
    way_pts = [tar_pos(:,1),cor,[190/3-x_dis,190/3-x_dis;30,110/3],tar_pos(:,3:4),[70-x_dis,64-x_dis;190/3,190/3],tar_pos(:,5:6)];
    
elseif (strcmp(type,'obs'))
    % obstacle way points
    c_set = inPara.c_set;
    r_set = inPara.r_set;
    theta_set = inPara.theta_set;
    rec_vert = inPara.rec_vert;
        
    obs_wp = cell(size(theta_set,1)+size(rec_vert,1),1);
    % generate way points for round obstacles
%     figure
    for ii = 1:size(theta_set,1)
        theta = theta_set{ii};
        cor = zeros(2,length(theta)); %coordinates
        for jj = 1:length(theta)
            cor(1,jj) = c_set(1,ii)+r_set(ii)*cos(theta(jj));
            cor(2,jj) = c_set(2,ii)+r_set(ii)*sin(theta(jj));
        end
        
%         hold on
%         fill(cor(1,:),cor(2,:),'r');

        obs_wp(ii) = {cor};
    end
    
    % generate rectangular obstacles
    for ii = 1:size(rec_vert,1)
       ll = rec_vert{ii}(:,1); % lower-left coordinate 
       ur = rec_vert{ii}(:,2); % upper-right coordinate 
       obs_wp(size(theta_set,1)+ii) = {[ll,[ur(1);ll(2)],ur,[ll(1);ur(2)]]};%
    end
    way_pts = obs_wp; 
end
end