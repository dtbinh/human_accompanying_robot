% Defines field class
% Field is assumed to be a rectange with endpoints=[xMin xMax yMin yMax]
classdef field
    
    properties
        targetNum; % specify target number
        targetPos; % specify target position, [x; y;]
        agentNum;
        endpoints; % endpoints of field, [xMin xMax yMin yMax]
    end
   
    methods
        function obj=field(endpoints,targetPos)
            if nargin > 0
                obj.endpoints = endpoints;
                obj.targetPos = targetPos;
                obj.targetNum = length(targetPos(1,:));
            end
        end
    end
end
        
        