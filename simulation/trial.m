load('../exp/SLAM_map/shared_room_map.mat')
occu_map = double(map.occu_map);
image(occu_map)
colormap(gray(256))