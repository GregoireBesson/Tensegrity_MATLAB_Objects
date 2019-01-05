% Function that returns random strings to actuate
%
% Inputs: - actuators (2 x 12) pair of strings indexes
%         - nbActuators to actuate
%         - rngParam for the repeatability
%         - actStringsIn (nbIndiv x nbActuation x 2*nbActuators),table 
%           containing the actuated strings so far
%         - indiv, individual counter
%         - actCounterIn, actuation counter
%
% Output: - actStringsOut (nbIndiv x nbActuation x 2*nbActuators), table
%           containing the updated strings to tension
%         - actCounterOut updated actuation counter

function [actStringsOut,actCounterOut] = randomStrings(actuators,nbActuators,rngParam,actStringsIn,indiv,actCounterIn)

% Draw nbActuators unique values from the integers 1 to 1.
randomIndexes = datasample(rngParam,1:12,nbActuators,'Replace',false); 

% vector containing which strings are going to be pulled
actuatedStrings = actuators(:,randomIndexes);
% reshape the 2 x nbAct matrix to a 1 x 2*nbAct Vector
actuatedStrings = reshape(actuatedStrings,[1,2*nbActuators]);

% increment the actuation counter
actCounterOut = actCounterIn + 1;

% add the new strings in the array
actStringsIn(indiv,actCounterOut,:) = actuatedStrings;

% update the output
actStringsOut = actStringsIn;


end