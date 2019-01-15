  function [actuatedSpringsTableOut,actCounterOut] = myDynamicsUpdate(tensStruct1, dynamicsPlot1, displayTimeInterval, actuatedSpringsIn, pretension, maxTension, l0, actuators,nbActuators,rngParam,indiv,actCounterIn,nbActuations,displaySimulation)
% This function will perform dynamics update each timestep.

%create some persistent variables for objects and structs
persistent tensStruct dynamicsPlot tspan i actuatedSpringsTableIn individual... 
           pretensionVector k lZero j m actuatorss nbAct rngPar actCounter...
           nbActMax displaySim

if nargin>1
    displaySim = displaySimulation;
    i = 0;          % counts the number of times this function is called
    j = 1;          % select wich couple of strings is actuated currently
    k = 1;          % incremented each time the pretension increases
    tensStruct = tensStruct1;
    dynamicsPlot = dynamicsPlot1;
    tspan = displayTimeInterval;
    actuatedSpringsTableIn = actuatedSpringsIn;
    actuatorss = actuators;
    nbAct = nbActuators;
    rngPar = rngParam;
    individual = indiv;
    actCounter = actCounterIn; 
    nbActMax = nbActuations;
    [~, m, ~] = size(actuatedSpringsTableIn);
    pretensionVector = [pretension:maxTension maxTension maxTension ...
             maxTension maxTension maxTension maxTension maxTension ... 
             maxTension maxTension maxTension maxTension maxTension ...
             pretension pretension pretension pretension pretension]; 
                                    %to hold the compressed position
    lZero = l0;
    
end

actuatedSpringsTableOut = actuatedSpringsTableIn;
actCounterOut = actCounter;

%%% Rest-length/tension controller %%%
i = i + 1;
if i > 50  % Start compression after a certain time.
    if (mod(i,2)==0) % increase tension every 20 loops --> jump in 50sec
        k = k + 1;
        % after reaching maxtension, release energy by resetting
        if (k > length(pretensionVector)) 
            % reshuffle random actuators
            if (actCounter <= nbActMax)
                if (strcmpi(tensStruct.selectionMode,'Random'))
                    [actuatedSpringsTableIn,actCounter] = randomStrings(actuatorss,nbAct,rngPar,actuatedSpringsTableIn,individual,actCounter);
                    actuatedSpringsTableOut = actuatedSpringsTableIn;
                    actCounterOut = actCounter;
                end
                % or go on the next couple of strings to be actuated
                if (m > 1)
                    j = j + 1;
                    if (j > m)
                        j = 1; % restart with first couple of strings when finished
                        actCounter = nbActMax + 1; %to stop the actuation
                    end
                end
                k = 1;
            end
        end
        if (actCounter <= nbActMax)
            % update only the actuated strings
            if (strcmpi(tensStruct.selectionMode,'Random'))
                tensStruct.simStruct.stringRestLengths(actuatedSpringsTableOut(individual,actCounter,:)) = ((100-pretensionVector(k))/100)*lZero;
            else
                tensStruct.simStruct.stringRestLengths(actuatedSpringsTableOut(individual,j,:)) = ((100-pretensionVector(k))/100)*lZero;
            end
        end
    end 
end

%%% End controller %%%

% Update nodes:
dynamicsUpdate(tensStruct, tspan);
dynamicsPlot.nodePoints = tensStruct.ySim(1:end/2,:);
if displaySim
updatePlot(dynamicsPlot);

drawnow  %plot it up
end

