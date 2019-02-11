% Author:   Gregoire Besson, EPFL student, gregoire.besson@alumni.epfl.ch
% Date:     Winter 2018-2019
%
% Intitializes the genome with random actuators
%
% Inputs: - nbActuationCycle
%         - AvlblActuators, list of available actuators
%         - nbUsedActuators, nb of actuators we want to actu simultaneously
%         - rngParam, for the repeatability        
%
% Output: - genesOut, (nbActuation x 12) genome output for nbAct of 1 indiv

function genesOut = initRandomGenome(nbActuationCycle,AvlblActuators,nbUsedActuators,rngParam)
    
    genes = zeros(nbActuationCycle,12);
    nbAvlAct = length(AvlblActuators);
    
    for c = 1:nbActuationCycle
        
        % Draw random indexes from the integers 1 to nbAvlAct.
        randomIndexes = datasample(rngParam,1:nbAvlAct,nbUsedActuators,'Replace',false);
        % indexes are used to find the actuators 
        usedActuators = AvlblActuators(randomIndexes);
        
        % set the motors to be actuated
        genes(c,usedActuators) = 1;
        
    end
    
    % update the output
    genesOut = genes;
    
end