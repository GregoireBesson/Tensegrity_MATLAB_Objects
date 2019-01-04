  function myDynamicsUpdate(tensStruct1, dynamicsPlot1, displayTimeInterval, actuatedSprings, pretension, maxTension, l0, actuators,nbActuators,rngParam)
% This function will perform dynamics update each timestep.

%create some persistent variables for objects and structs
persistent tensStruct dynamicsPlot tspan i actuatedSpringsVec ... 
           pretensionVector k lZero j m actuatorss nbAct rngPar

if nargin>1
    i = 0;          % counts the number of times this function is called
    j = 1;          % select wich couple of strings is actuated currently
    k = 1;          % incremented each time the pretension increases
    tensStruct = tensStruct1;
    dynamicsPlot = dynamicsPlot1;
    tspan = displayTimeInterval;
    actuatedSpringsVec = actuatedSprings;
    actuatorss = actuators;
    nbAct = nbActuators;
    rngPar = rngParam;
    [m, ~] = size(actuatedSpringsVec);
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
    if (mod(i,2)==0) % increase tension every 20 loops --> jump in 50sec
        k = k + 1;
        % after reaching maxtension, release energy by resetting
        if (k > length(pretensionVector)) 
            % reshuffle random actuators
            actuatedSpringsVec = randomStrings(actuatorss,nbAct,rngPar);
            % or go on the next couple of strings to be actuated
            if (m > 1)
                j = j + 1;
                % restart with first couple of strings when finished
                if (j > m)
                   j = 1; 
                end
            end
            k = 1;
        end
    %newRestLengths =  tensStruct.simStruct.stringRestLengths;
    %newRestLengths(ActuatedStrings) = ((100-pretensionVector(k))/100)*norm(tensStruct.nodePoints(1,:)-tensStruct.nodePoints(7,:));
    
    % update only the actuated strings
    tensStruct.simStruct.stringRestLengths(actuatedSpringsVec(j,:)) = ((100-pretensionVector(k))/100)*lZero;
    end 
end

%%% End controller %%%

% Update nodes:
dynamicsUpdate(tensStruct, tspan);
dynamicsPlot.nodePoints = tensStruct.ySim(1:end/2,:);
updatePlot(dynamicsPlot);

drawnow  %plot it up
end

