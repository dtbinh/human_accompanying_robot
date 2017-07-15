function [action,actionIndex]=getActionFromInput(agent,field)
    inputAction=input('Choose human action: up, down, left, right\n','s');
    x = agent.currentPos(1);
    y = agent.currentPos(2);
    psi = agent.currentPos(3);
    xMin = field.endpoints(1);
    xMax = field.endpoints(2);
    yMin = field.endpoints(3);
    yMax = field.endpoints(4);
    tolerance = 1e-4;
    legalActions = returnLegalActions(agent,field);
    
    switch lower (inputAction)
        case 'up'
            newX=x;
            newY=y+agent.maxV;
            newPsi=pi/2;
            
        case 'down'
            newX=x;
            newY=y-agent.maxV;
            newPsi=3/2*pi;
            
        case 'left'
            newX=x-agent.maxV;
            newY=y;
            newPsi=pi;
            
        case 'right'
            newX=x+agent.maxV;
            newY=y;
            newPsi=0;
        
        otherwise
            display('Invalid Action, please input again\n');
            [action,actionIndex] = getActionFromInput(agent,field);
            return
    end
    
    newX = min(newX,xMax);
    newX = max(newX,xMin);
    newY = min(newY,yMax);
    newY = max(newY,yMin);
    action = [newX-x;newY-y;newPsi-psi];
%     [~,actionIndex] = ismember(action',legalActions','rows');
    for actionIndex = 1:length(legalActions(1,:))
        if abs(action-legalActions(:,actionIndex))<tolerance
            break
        elseif actionIndex == length(legalActions(1,:))
            error('action not found in legal actions')
        end
    end
end
        