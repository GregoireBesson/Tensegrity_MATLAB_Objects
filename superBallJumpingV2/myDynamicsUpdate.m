% Author:     ???
% ModifiedBy: Gregoire Besson, EPFL student, gregoire.besson@alumni.epfl.ch
% Date:       Winter 2018-2019
%
% This function will perform dynamics update each timestep.

function myDynamicsUpdate(tensStruct1, dynamicsPlot1, displayTimeInterval, pretension, maxTension, l0, actuators, indiv,nbActuationCycle,displaySimulation, genesIn,g,firstStringsToActuate)

%create some persistent variables for objects and structs
persistent tensStruct dynamicsPlot tspan i individual ... 
           pretensionVector k lZero actuatorss stringsToActuate...
           nbActMax displaySim genes actCycleCounter genNumber

if nargin>1
    displaySim = displaySimulation;
    i = 0;          % counts the number of times this function is called
    k = 1;          % incremented each time the pretension increases
    actCycleCounter = 1;
    tensStruct = tensStruct1;
    dynamicsPlot = dynamicsPlot1;
    tspan = displayTimeInterval;
    actuatorss = actuators;
    individual = indiv;
    genes = genesIn;
    nbActMax = nbActuationCycle;
    pretensionVector = [pretension:maxTension maxTension maxTension ...
             maxTension maxTension maxTension maxTension maxTension ... 
             maxTension maxTension maxTension maxTension maxTension ...
             pretension pretension pretension pretension pretension]; 
                                    %to hold the compressed position
    lZero = l0;
    genNumber = g;
    stringsToActuate = firstStringsToActuate;
    
end

%%% Rest-length/tension controller %%%

i = i + 1;
% Start compression after a certain time.
if i > 50  
    % increase tension every 20 loops --> jump in 50sec
    if (mod(i,2)==0)   
            k = k + 1;  
            % after reaching maxtension, release energy by resetting
            if ( k > length(pretensionVector) )
                % reset the tension vector counter
                k = 1;
                % increment the Actuation Cycle
                actCycleCounter = actCycleCounter + 1;
                % safety to stop the actuation cycle
                if (actCycleCounter > nbActMax)
                    actCycleCounter = 1;
                    i = -10000;
                end
                % update the strings to actuate
                stringsToActuate = genes2strings(genes,genNumber,individual,actCycleCounter,actuatorss);
            end
            % update tension of the strings to actuate
            tensStruct.simStruct.stringRestLengths(stringsToActuate) = ((100-pretensionVector(k))/100)*lZero;      
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

