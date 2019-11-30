% Author:     Jeff Friesen
% ModifiedBy: Gregoire Besson, EPFL student, gregoire.besson@alumni.epfl.ch
% Date:       Winter 2018-2019
%
% This function will perform dynamics update each timestep.

function myDynamicsUpdate(tensStruct1, dynamicsPlot1, displayTimeInterval, pretension, maxTension, l0, actuators, indiv,nbActuationCycle,displaySimulation, genesIn,g,firstStringsToActuate)

%create some persistent variables for objects and structs
persistent tensStruct dynamicsPlot tspan i individual tensionVector k l ...
           lZeroActive actuatorss stringsToActuate nbActMax displaySim genes ...
           actCycleCounter genNumber tensionVectorNoRelease stringsToHold ...
           holdCounter

if nargin>1
    displaySim = displaySimulation;
    i = 0;          % counts the number of times this function is called
    k = 1;          % incremented each time the pretension increases
    l = 1;          % same but used for the hold maxtension 
    actCycleCounter = 1;
    tensStruct = tensStruct1;
    dynamicsPlot = dynamicsPlot1;
    tspan = displayTimeInterval;
    actuatorss = actuators;
    individual = indiv;
    genes = genesIn;
    nbActMax = nbActuationCycle;
    % loockup table for increasing tension then release
    tensionVector = [pretension:maxTension maxTension maxTension ...
             maxTension maxTension maxTension maxTension maxTension ... 
             maxTension maxTension maxTension maxTension maxTension ...
             pretension pretension pretension pretension pretension ]; 
    % loockup table for increasing and holding the max tension                               
    tensionVectorNoRelease = [pretension:maxTension maxTension maxTension ...
             maxTension maxTension maxTension maxTension maxTension ... 
             maxTension maxTension maxTension maxTension maxTension ...
             maxTension maxTension maxTension maxTension maxTension];                                
    lZeroActive = l0;
    genNumber = g;
    stringsToActuate = firstStringsToActuate;
    %stringsToHold = [9 10 ; nice locomotion result
%                      5  6 ;
%                      9 10 ;
%                      5  6];
%                    loc1   loc1 loc2 loc3
    stringsToHold = [ 15 ;  %15   %15  %24
                      6 ;   %6 ; %4;  %17
                      15 ]; %15   %15  %24
    holdCounter = 1;
    
end

%%% Rest-length/tension controller %%%

i = i + 1;
% Start compression after a certain time.
if i > 50  
    % increase tension every 3 loops --> jump in 50sec
    if (mod(i,3)==0)   
            k = k + 1; 
            l = l + 1;
            % after reaching maxtension, release energy by resetting
            if ( k > length(tensionVector) )
                % reset the tension vector counter
                k = 1;
                % but not of the hold tension counter
                l = 1;
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
%             if ( l > length(tensionVectorNoRelease) )
%                 % reset the hold tension vector counter and release
%                 l = 1;
%                 holdCounter = holdCounter + 1;
%             end
            % update tension of the strings to actuate
            tensStruct.simStruct.stringRestLengths(stringsToActuate) = ((100-tensionVector(k))/100)*lZeroActive;
            % and update tension of the specific strings to hold
            %tensStruct.simStruct.stringRestLengths(stringsToHold(actCycleCounter,:)) = ((100-tensionVectorNoRelease(l))/100)*lZeroActive;

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

