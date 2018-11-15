  function myDynamicsUpdate(tensStruct1, dynamicsPlot1, displayTimeInterval, actuatedSprings, pretension, maxTension, l0)
% This function will perform dynamics update each timestep.

%create some persistent variables for objects and structs
persistent tensStruct dynamicsPlot tspan i actuatedSpringsVec pretensionVector k lZero

if nargin>1
    i = 0;          % counts the number of times this function is called
    k = 1;          % incremented each time the pretension increases
    tensStruct = tensStruct1;
    dynamicsPlot = dynamicsPlot1;
    tspan = displayTimeInterval;
    actuatedSpringsVec = actuatedSprings;
    pretensionVector = [pretension:maxTension maxTension maxTension ...
             maxTension maxTension maxTension maxTension maxTension ... 
             maxTension maxTension maxTension maxTension maxTension ...
             pretension pretension pretension pretension pretension]; 
                                    %to hold the compressed position
    lZero = l0;
end

%%% Optional rest-length controller %%%
i = i + 1;
if i > 50  % Start compression after a certain time.
    if (mod(i,3)==0) % increase tension every 2 loops
        k = k + 1;
        % after reaching maxtension, release energy by resetting
        if (k > length(pretensionVector)) 
            k = 1;
        end
    %newRestLengths =  tensStruct.simStruct.stringRestLengths;
    %newRestLengths(ActuatedStrings) = ((100-pretensionVector(k))/100)*norm(tensStruct.nodePoints(1,:)-tensStruct.nodePoints(7,:));
    
    % update only the actuated strings
    tensStruct.simStruct.stringRestLengths(actuatedSpringsVec) = ((100-pretensionVector(k))/100)*lZero;
    end 
end

%%% End controller %%%

% Update nodes:
dynamicsUpdate(tensStruct, tspan);
dynamicsPlot.nodePoints = tensStruct.ySim(1:end/2,:);
updatePlot(dynamicsPlot);

drawnow  %plot it up
end

