% Function that returns random strings to actuate
%
% Inputs: - actuators, (2 x 12) pair of strings indexes
%         - nbActuators, to actuate
%         - rngParam, for the repeatability
%         - actStringsIn, (nbIndiv x nbActuation x 2*nbActuators),table 
%           containing the actuated strings so far
%         - genesIn, (nbGen x nbIndiv x nbActuation x 12) genome input
%         - indiv, individual counter
%         - g, generation counter
%         - actCounterIn, actuation counter
%
% Output: - actStringsOut, (nbIndiv x nbActuation x 2*nbActuators), table
%           containing the updated strings to tension
%         - genesOut, (nbGen x nbIndiv x nbActuation x 12) genome output
%         - actCounterOut, updated actuation counter

function [actStringsOut,genesOut,actCounterOut] = randomStrings(actuators,nbActuators,rngParam,actStringsIn,genesIn,indiv,g,actCounterIn)

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