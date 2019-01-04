% Function that returns random strings to actuate
%
% Inputs: - actuatedPairStrings (2 x 12) pair of strings indexes
%         - nbActuators to actuate
%         - stringNumber for the repeatability
%
% Output: - actuatedStrings (1 x nbActuators*2) containing the strings to
%           tension

function actuatedStrings = randomStrings(actuatedPairStrings,nbActuators,stringNumber)

%Draw nbActuators unique values from the integers 1 to 1.
randomIndexes = datasample(stringNumber,1:12,nbActuators,'Replace',false); 

% vector containing which strings are going to be pulled
actuatedStrings = zeros(2,nbActuators);

for i = 1:length(randomIndexes)
    actuatedStrings(:,i) = actuatedPairStrings(:,randomIndexes(i));
end

actuatedStrings = reshape(actuatedStrings,[1,2*nbActuators]);

end