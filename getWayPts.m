function way_pts = getWayPts(tar_pos,scale,type)
if (strcmp(type,'h'))
    % human way points
    % generate way points for human trajectory around the round obstacle
    center = [220;40]; %center of circle
    r = 15; %radius
    theta = -3/4*pi:pi/8:pi/4;
    cor = zeros(2,length(theta)); %coordinates
    for ii = 1:length(theta)
        cor(1,ii) = center(1)+r*cos(theta(ii));
        cor(2,ii) = center(2)+r*sin(theta(ii));
    end
    figure
    hold on
    fill(cor(1,:),cor(2,:),'r');
    way_pts = [tar_pos(:,1),cor*scale,[190,190;90,110]*scale,tar_pos(:,3:4),[210,190;190,190]*scale,tar_pos(:,5)];
    
elseif (strcmp(type,'obs'))
    % obstacle way points
    % generate way points for the round obstacle
    center = [220;40]; %center of circle
    r = 10; %radius
    theta = 0:pi/8:2*pi;
    cor = zeros(2,length(theta)); %coordinates
    for ii = 1:length(theta)
        cor(1,ii) = center(1)+r*cos(theta(ii));
        cor(2,ii) = center(2)+r*sin(theta(ii));
    end
%     figure
%     hold on
%     fill(cor(1,:),cor(2,:),'r');
    way_pts = cor*scale;
end
end