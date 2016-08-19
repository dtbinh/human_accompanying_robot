% load('../exp/SLAM_map/shared_room_map.mat')
% occu_map = double(map.occu_map);
% image(occu_map)
% colormap(gray(256))

max_num = 20;
len = size(traj,2);
t = floor(len/max_num);
traj_p = [];

for ii = 1:t
    if ii < t
        tmp_pts = traj(:,max_num*(ii-1)+1:max_num*ii);
    else
        tmp_pts = traj(:,max_num*(ii-1)+1:end);
    end
    % polynomial fitting
%     figure
%     plot(tmp_pts(1,:),tmp_pts(2,:),'o')
%     hold on
    p1 = fit(tmp_pts(1,:)',tmp_pts(2,:)','pchipinterp');
%     plot(p1,'r')
    y1 = p1(tmp_pts(1,:));
    traj_p = [traj_p,[tmp_pts(1,:);y1']];
%     hold on
%     plot(tmp_pts(1,:),y1','g')
%     p2 = spline(tmp_pts(1,:)',tmp_pts(2,:)');
%     y2 = ppval(p2,tmp_pts(1,:)');
%     plot(tmp_pts(1,:)',y2,'r')
end
figure
plot(traj(1,:),traj(2,:),'o')
hold on
plot(traj_p(1,:),traj_p(2,:),'g')