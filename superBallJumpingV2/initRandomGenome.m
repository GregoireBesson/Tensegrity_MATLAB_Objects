% Author:   Gregoire Besson, EPFL student, gregoire.besson@alumni.epfl.ch
% Date:     Winter 2018-2019
%
% Intitializes the genome with random actuators
%
% Inputs: - nbActuationCycle
%         - nbActuators, to actuate
%         - nbIndiv, size of poupulation
%         - rngParam, for the repeatability        
%         - indiv, individual counter
%
% Output: - genesOut, (nbActuation x 12) genome output for nbAct of 1 indiv

function genesOut = initRandomGenome(nbActuationCycle,nbActuators,rngParam)
    
    genes = zeros(nbActuationCycle,12);

    for c = 1:nbActuationCycle
        
        % Draw nbActuators unique values from the integers 1 to 1.
        randomIndexes = datasample(rngParam,1:12,nbActuators,'Replace',false);
        
        % set the motors to be actuated
        genes(c,randomIndexes) = 1;
        
    end
    
    % update the output
    genesOut = genes;
    
end