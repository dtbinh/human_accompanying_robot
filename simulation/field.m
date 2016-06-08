% Defines field class
% Field is assumed to be a rectange with endpoints=[xMin xMax yMin yMax]
classdef field
    
    properties
        tar_num; % specify target number
        tar_pos; % specify target position, [x; y;]
        agnt_num;
        end_pts; % end_pts of field, [xMin xMax yMin yMax]
        obs_info; % information of obstacles
    end
   
    methods
        function this=field(end_pts,tar_pos)
            if nargin > 0
                this.end_pts = end_pts;
                this.tar_pos = tar_pos;
                this.tar_num = length(tar_pos(1,:));
            end
        end
        
        function this = processOccuMap(this,map)
            % process the raw occupancy map by fitting polygons to obstacles
            % cluster
            clst = kmeans(map.)
            % fit convex hull
            
        end
    end
end
        
        