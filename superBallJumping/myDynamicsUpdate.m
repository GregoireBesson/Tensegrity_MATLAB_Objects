function myDynamicsUpdate(tensStruct1, dynamicsPlot1, displayTimeInterval)
% This function will perform dynamics update each timestep.

%create some persistent variables for objects and structs
persistent tensStruct dynamicsPlot tspan

if nargin>1
    tensStruct = tensStruct1;
    dynamicsPlot = dynamicsPlot1;
    tspan = displayTimeInterval;
end

% Update nodes:
dynamicsUpdate(tensStruct, tspan);
dynamicsPlot.nodePoints = tensStruct.ySim(1:end/2,:);
updatePlot(dynamicsPlot);

drawnow  %plot it up
end