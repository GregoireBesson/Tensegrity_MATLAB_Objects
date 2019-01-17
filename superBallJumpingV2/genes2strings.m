function stringsToActuate = genes2strings(genes,genNumber,individual,actCycleCounter,actuators,nbActuators)
    
    % find wich actuators has a positive gene
    motorIndexes =  genes(genNumber,individual,actCycleCounter,:) == 1 ;
    % deduce the corresponding strings to actuate
    stringIndexes = actuators(:,motorIndexes);
    % reshape to the output format
    stringsToActuate = reshape(stringIndexes,1,2*nbActuators);

end